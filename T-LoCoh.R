#安装必要的包----------
# 更换为合适的镜像
options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))
library(sp)
library(raster)
library(geosphere)
library(sf)
library(move)
library(FNN)
library(png)
library(XML)
# 从 R - Forge 下载 tlocoh 包的源码
download.file("http://R-Forge.R-project.org/src/contrib/tlocoh_1.40.07.tar.gz", "tlocoh.tar.gz")
#为什么tlocoh每次要重新下
# 从源码安装包
install.packages("tlocoh.tar.gz", repos = NULL, type = "source")
library(tlocoh)
library(dplyr)
library(base)
library(stats)
library(rgdal)

setwd("C:/data/R/github")
T1data <- read.csv("T1/T1data.csv")
# 查看数据
class(T1data)
head(T1data)
plot(T1data[, c("location.long","location.lat")], pch=20)

#转换成投影坐标
T1data.sp.latlong <- SpatialPoints(T1data[, c("location.long","location.lat")], proj4string=CRS("+proj=longlat +ellps=WGS84")) 
T1data.sp.utm <- spTransform(T1data.sp.latlong, CRS("+proj=utm +north +zone=52 +ellps=WGS84")) 
T1data.mat.utm <- coordinates(T1data.sp.utm) 
colnames(T1data.mat.utm) <- c("x","y")
head(T1data.mat.utm)

# 时间转换
class(T1data$t)
T1data.gmt <- as.POSIXct(T1data$t, tz="UTC") 
local.tz <- "Asia/Shanghai" 
T1data.localtime <- as.POSIXct(format(T1data.gmt, tz=local.tz), tz=local.tz) 
T1data.localtime[1:3]

#创建 locoh - xy 对象
#将动物运动轨迹的位置、时间戳、个体标识等数据整合在一起，便于统一管理和操作,尤其是多个体的数据
T1data.lxy <- xyt.lxy(xy=T1data.mat.utm, dt=T1data.localtime, id="T1") 

#数据清洗和查看
##数据基本信息，整体特征。查看重复数值（相同或不同时间戳和位置信息）
summary(T1data.lxy)
##检查位置。不同颜色即不同时间段的运动
plot(T1data.lxy)
##直方图，步长
hist(T1data.lxy)
##检查采样间隔（是否一致）
##检查按日期划分的采样频率图
lxy.plot.freq(T1data.lxy, deltat.by.date=T)
lxy.plot.freq(T1data.lxy, cp=T)
##将短于中位采样间隔 20% 的点视为 “爆发点” 并稀释，若数据中存在密集采样段（如GPS故障），需适当降低阈值（如thresh=0.1）；若数据均匀，可提高阈值。
T1data.lxy <- lxy.thin.bursts(T1data.lxy, thresh=0.2)

#子集操作
##数据集大且无需构建全部hulls时，比如仅试图创建利用率分布。
##抽样间隔=3，用于减少计算量。信息损失风险，特别是信息稀疏，行为复杂
T1data.lxy.every3rd <- lxy.subset(T1data.lxy, idx = seq(from=1, to=5775, by=3)) 
summary(T1data.lxy.every3rd) 

#确定时间平衡
T1data.lxy <- lxy.ptsh.add(T1data.lxy)
lxy.plot.pt2ctr(T1data.lxy)
##确定s值（控制时间间隔和空间间隔在临近点选择的权重）
###s=0，仅空间距离决定邻近点，时间无影响（传统LoCoH），低s或=0.0005
###s＞0，越大，时间主导权越大。高s尺度：0.01，中s常为0.003（即时空平衡）
lxy.plot.sfinder(T1data.lxy) 
##指定要分析的时间间隔，t=s
lxy.plot.sfinder(T1data.lxy, delta.t=3600*c(12,24,36,48,54,60))

#识别临近点
##k=25表示每个点选择25个最近邻构建hull。
T1data.lxy <- lxy.nn.add(T1data.lxy, s=0.003, k=25)
##使用s值的范围构建hulls
T1data.lxy <- lxy.nn.add(T1data.lxy, s=c(0.0003, 0.003, 0.03, 0.3), k=25) 
##评估不同k值的等值线面积和边缘比。
##k值较大：hull覆盖范围广，可能包含未使用区域（Type I错误）
##k值较小：hull碎片化，可能遗漏核心区域（Type II错误）。
lxy.plot.mtdr(T1data.lxy, k=10) 
lxy.plot.tspan(T1data.lxy, k=10)

#保存数据
lxy.save(T1data.lxy, dir=".") 

#创建 hullsets
##尝试k值为6,9,12,...,24的多种组合,选择最优值
T1data.lhs <- lxy.lhs(T1data.lxy, k=3*2:8, s=0.003)
summary(T1data.lhs, compact=T) 
#构建等值线
T1data.lhs <- lhs.iso.add(T1data.lhs)
plot(T1data.lhs, iso=T, record=T, ufipt=F)
plot(T1data.lhs, iso=T, k=15, allpts=T, cex.allpts=0.1, col.allpts="gray30", ufipt=F)
lhs.plot.isoarea(T1data.lhs)
lhs.plot.isoear(T1data.lhs)
T1data.lhs.k15 <- lhs.select(T1data.lhs, k=15)

#a - 方法操作
T1data.lxy <- lxy.nn.add(T1data.lxy, s=0.003, a=auto.a(nnn=15, ptp=0.98)) 
T1data.lxy <- lxy.nn.add(T1data.lxy, s=0.003, a=15000) 
T1data.lhs.amixed <- lxy.lhs(T1data.lxy, s=0.003, a=4:15*1000, iso.add=T)
lhs.plot.isoarea(T1data.lhs.amixed) 
lhs.plot.isoear(T1data.lhs.amixed)

#计算额外的 hull 指标
T1data.lhs.k15 <- lhs.ellipses.add(T1data.lhs.k15) 
summary(T1data.lhs.k15) 
plot(T1data.lhs.k15, hulls=T, ellipses=T, allpts=T, nn=T, ptid="auto")
T1data.lhs.k15 <- lhs.visit.add(T1data.lhs.k15, ivg=3600*12) 
summary(T1data.lhs.k15) 

#查看hull指标
lhs.iso.add(T1data.lhs.k15, sort.metric="ecc") 
plot(T1data.lhs.k15, iso=T, iso.sort.metric="ecc")
hist(T1data.lhs.k15, metric="nsv") 
plot(T1data.lhs.k15, hpp=T, hpp.classify="nsv", ivg=3600*12, col.ramp="rainbow")
