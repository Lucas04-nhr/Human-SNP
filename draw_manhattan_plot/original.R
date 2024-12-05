## ----pressure1, echo=FALSE-------------------------------------------------------------
library(tidyverse)
library(ggforce)
library(ggprism)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)
df_all<-read.table("result.P16.assoc.linear.adjusted",col.names = c("CHR","SNP","UNADJ","GC","BONF","HOLM","SIDAK_SS","SIDAK_SD","FDR_BH","FDR_BY","Bacteria"),skip = 1)#导入文件
df_all$BP <- sapply(df_all$SNP, function(x) {
  parts <- strsplit(x, ":")
  return(parts[[1]][2])  # 返回分割后的第二部分
})



## ----pressure2, echo=FALSE-------------------------------------------------------------

sorted_df_all <- df_all %>%
  arrange(CHR, BP)
sorted_df<-sorted_df_all[,c("CHR","BP","BONF","Bacteria")]


## ----pressure3, echo=FALSE-------------------------------------------------------------
sorted_df$BP<-as.numeric(sorted_df$BP)
name<-unique(sorted_df_all$CHR)
sorted_df_all$CHR<-factor(sorted_df_all$CHR, levels =c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","26",name[25:length(name)]))


## ----pressure4, echo=FALSE-------------------------------------------------------------


lengths <- sorted_df %>%
  group_by(CHR) %>%
  summarise(
    min_position = min(BP),
    max_position = max(BP),
    length_approx = max_position - min_position + 1 # 加1是因为位置是离散的，包含两端点
  )



## ----pressure5, echo=FALSE-------------------------------------------------------------

lengths_chr <- lengths %>% filter(CHR %in% c(1:23,26))
lengths_chr$CHR<-as.numeric(lengths_chr$CHR)
lengths_chr <- lengths_chr %>% arrange(CHR)


## ----pressure6, echo=FALSE-------------------------------------------------------------

lengths_other <- lengths[!(lengths$CHR %in% c(1:23,26)),]
lengths_other <- lengths_other %>% arrange(CHR)


## ----pressure7, echo=FALSE-------------------------------------------------------------

len<-rbind(lengths_chr,lengths_other)
len<-len %>% mutate(Sum=cumsum(length_approx)-length_approx)


## ----pressure8, echo=FALSE-------------------------------------------------------------

sorted_df_len <- full_join(sorted_df, len[,c("CHR","Sum")], by = "CHR")
df <- sorted_df_len %>% mutate(.,x=as.numeric(BP)+Sum)
df <- df%>%
  mutate(.,logP=-log10(BONF))


## ----pressure9, echo=FALSE-------------------------------------------------------------

# 将CHR中的23替换为X, 26替换为Y, 其他替换为Others
old_elements <- c("23", "26")
new_elements <- c("X", "Y")

df$CHR <- ifelse(df$CHR %in% old_elements, new_elements, df$CHR)


old_elements <- name[25:length(name)]
new_element <- "Other"
 
# 使用基本的R语法进行替换
df$CHR <- ifelse(df$CHR %in% old_elements, new_element, df$CHR)



## ----pressure10, echo=FALSE------------------------------------------------------------
chromosome_colors <- c(
  rep(c('#1F77B4','#FF7F0C','#2BA02B','#D62628','#9467BD','#8C564B','#7F7F7F','#E477C2','#BDBD21', '#17BECF'),3))
df_CHR<-df %>% filter(CHR%in% c(1:23,26))
v0<-sort(unique(df_CHR$Sum))
v0<-c(v0,0)
v0[length(v0)]<-(v0[length(v0)-1]+as.numeric(lengths_chr[nrow(lengths_chr),4]))
v1<-v0
for (i in 1:length(v1)) {
  v1[i]<-(as.numeric(v0[i])+as.numeric(v0[i+1]))/2}


## ----pressure11, echo=FALSE------------------------------------------------------------
a<-len[nrow(len),"length_approx"]
a<-as.numeric(a)
v1[length(v1)]<-(v1[length(v1)-1]+max(len$Sum)+a)/2


## ----pressure12, echo=FALSE------------------------------------------------------------
df$Bacteria<-as.character(df$Bacteria)
df <- df %>%
  group_by(Bacteria) %>%
  mutate(rank_p = rank(BONF, ties.method = "first")) %>% # 对P进行排名
  mutate(Bacteria_new = if_else(rank_p == 1, Bacteria, NA_character_)) %>% # 保留P最小的行的Bacteria值，其他变成NA
  ungroup() %>%
  select(-rank_p)



## ----pressure13, echo=FALSE------------------------------------------------------------
#修改CHR的levels
df$CHR<-factor(df$CHR, levels =c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","26","Other"))


## ----pressure14, echo=FALSE------------------------------------------------------------

p<-ggplot(df, aes(x = x, y = logP, color = as.factor(CHR))) +
  geom_point(alpha = 0.7, size = 1.5) +  # 绘制点图
  scale_color_manual(values = chromosome_colors) +  # 自定义颜色
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(df$logP) + 1)) +  # 调整Y轴范围，确保没有间隙
  theme_minimal() +  # 使用简洁的主题
  theme(legend.position = "none",  # 隐藏图例
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),  # 设置轴线粗细（单位为pt，即点）
  axis.ticks = element_line(size = 0.5),
  axis.ticks.length.y = unit(0.3, "cm") ,
  axis.ticks.length.x = unit(0, "cm")) + 
  geom_hline(yintercept = c(5, 6), color = c('blue', 'red'), linetype = c('dashed', 'dotted')) +  # 添加阈值线
  ggtitle("Manhattan Plot")+  # 添加标题
  labs(x="CHR",y = "-log10(P-value)", 
       title = "Manhattan Plot")+
    geom_text_repel(aes(label = df$Bacteria_new), color = "black", size = 3,force = 20, # 增加排斥力
    point.padding = 5)+scale_x_continuous(breaks =v1,labels =c(1:22,"X","Y","Other"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p


## ----pressure15, echo=FALSE------------------------------------------------------------
ggsave("p_16.pdf", plot = p, width = 10, height =4, units = "in", dpi = 300)

