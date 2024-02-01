library(readr)
library(scales)

diet = read.delim('../diet_aav_pc_concentrations.csv', stringsAsFactors = F, sep = ',')

for (i in 1:nrow(diet)){
  diet$av.wt[i] =  mean(as.numeric(diet[i,grep('WT',colnames(diet))]))
  diet$av.ko1[i] =  mean(as.numeric(diet[i,grep('KOS\\(1\\)',colnames(diet))]))
  diet$av.ko2[i] =  mean(as.numeric(diet[i,grep('KOS\\(2\\)',colnames(diet))]))
}

diet$length = sapply(strsplit(rownames(diet), ' '), '[', 3)                       
diet$dbonds = sapply(strsplit(diet$length, ':'), '[', 2)                       
diet$length = gsub('C', '', sapply(strsplit(diet$length, ':'), '[', 1))                      

diet$type = sapply(strsplit(rownames(diet), ' '), '[', 2)                       

colnames(diet)
diet.aa = diet[grep('aa', rownames(diet)), c(17:22)]
diet.ae = diet[grep('ae', rownames(diet)), c(17:22)]

### INDIVIDUAL PLOT - PC_aa class
pdf('Distribution_of_PC_aa_class_DIET_AAV_INDIVID.pdf', w = 14, h = 5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0), mfrow=c(1,1)) 
par(mar = c(4,4,4,7), mgp = c(3,1,0), mfrow=c(2,3)) 

plot(diet.aa$length, diet.aa$dbonds, cex = diet.aa$av.wt*0.014+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'WT', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('A')), bty = "n", inset = c(-0.2, -0.2), xpd = TRUE)

axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','50','100','  200'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.05,50,100,200)*0.014+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
# dev.off()


# pdf('Distribution_of_PC_aa_class_DIET_AAV_INDIVID_KO.pdf', w = 6, h = 3.5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0), mfrow=c(1,1)) 
plot(diet.aa$length, diet.aa$dbonds, cex = diet.aa$av.ko1*0.014+0.6, col = '#663333', bg = alpha('#663333', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'KOS(1)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('B')), bty = "n", inset = c(-0.2, -0.2), xpd = TRUE)
axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','50','100','  200'),col='#663333',pt.bg = alpha('#663333', 0.5),
       pch=21,pt.cex=c(0.05,50,100,200)*0.014+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
# dev.off()


# pdf('Distribution_of_PC_aa_class_DIET_AAV_INDIVID_KOS.pdf', w = 6, h = 3.5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0), mfrow=c(1,1)) 
plot(diet.aa$length, diet.aa$dbonds, cex = diet.aa$av.ko2*0.014+0.6, col = '#065c69', bg = alpha('#065c69', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'KOS(2)', bty = "n", inset = c(-0.02, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('C')), bty = "n", inset = c(-0.2, -0.2), xpd = TRUE)
axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','50','100','  200'),col='#065c69',pt.bg = alpha('#065c69', 0.5),
       pch=21,pt.cex=c(0.05,50,100,200)*0.014+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
# dev.off()


par(mar = c(6,4,4,1))
b = barplot(diet.aa$av.wt/sum(diet.aa$av.wt)*100, col = alpha('#b4b4b4',0.8), yaxt = 'n', ylab = 'Percentage')
legend('topleft', legend = expression(bold('D')), bty = "n", inset = c(-0.15, -0.25), xpd = TRUE)
axis(2, las = 2)
axis(1, b, labels = rownames(diet.aa), las = 2)
legend('topleft', legend = 'WT', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)

b = barplot(diet.aa$av.ko1/sum(diet.aa$av.ko1)*100, col = alpha('#663333', 0.8), yaxt = 'n', ylab = 'Percentage')
legend('topleft', legend = expression(bold('E')), bty = "n", inset = c(-0.15, -0.25), xpd = TRUE)
axis(2, las = 2)
axis(1, b, labels = rownames(diet.aa), las = 2)
legend('topleft', legend = 'KOS(1)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)

b = barplot(diet.aa$av.ko2/sum(diet.aa$av.ko2)*100, col = alpha('#065c69', 0.8), yaxt = 'n', ylab = 'Percentage')
legend('topleft', legend = expression(bold('F')), bty = "n", inset = c(-0.15, -0.25), xpd = TRUE)
axis(2, las = 2)
axis(1, b, labels = rownames(diet.aa), las = 2)
legend('topleft', legend = 'KOS(2)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)

dev.off()


### OVERLAY PLOT - PC_aa class
pdf('Distribution_of_PC_aa_class_DIET_AAV_OVERLAY.pdf', w = 12, h = 2.5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0),mfrow=c(1,1)) 
par(mar = c(4,4,4,7), mgp = c(3,1,0),mfrow=c(1,3)) 

plot(diet.aa$length, diet.aa$dbonds, cex = diet.aa$av.wt*0.014+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'WT + KOS(1)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('A')), bty = "n", inset = c(-0.25, -0.2), xpd = TRUE)

axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','50','100','  200'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.05,50,100,200)*0.014+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
legend('bottomright', legend = paste0(c('WT', 'KOS(1)')),col=c('#b4b4b4', '#663333'),
       pt.bg = c(alpha('#b4b4b4', 0.5), alpha('#663333', 0.5)),pch=21,pt.cex=2.75,
       cex = 0.8,title = '', xpd = TRUE, inset = c(-0.19,0), bty = 'n')
# dev.off()
points(diet.aa$length, diet.aa$dbonds, cex = diet.aa$av.ko1*0.014+0.6, col = '#663333', bg = alpha('#663333', 0.5), pch = 21)

plot(diet.aa$length, diet.aa$dbonds, cex = diet.aa$av.wt*0.014+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'WT + KOS(2)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('B')), bty = "n", inset = c(-0.25, -0.2), xpd = TRUE)

axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','50','100','  200'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.05,50,100,200)*0.014+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
legend('bottomright', legend = paste0(c('WT', 'KOS(2)')),col=c('#b4b4b4', '#065c69'),
       pt.bg = c(alpha('#b4b4b4', 0.5), alpha('#065c69', 0.5)),pch=21,pt.cex=2.75,
       cex = 0.8,title = '', xpd = TRUE, inset = c(-0.21,0), bty = 'n')
# dev.off()
points(diet.aa$length, diet.aa$dbonds, cex = diet.aa$av.ko2*0.014+0.6, col = '#065c69', bg = alpha('#065c69', 0.5), pch = 21)


plot(diet.aa$length, diet.aa$dbonds, cex = diet.aa$av.wt*0.014+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'WT + KOS(1) + KOS(2)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('C')), bty = "n", inset = c(-0.25, -0.2), xpd = TRUE)

axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','50','100','  200'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.05,50,100,200)*0.014+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
legend('bottomright', legend = paste0(c('WT', 'KOS(1)', 'KOS(2)')),col=c('#b4b4b4', '#663333','#065c69'),
       pt.bg = c(alpha('#b4b4b4', 0.5), alpha('#663333', 0.5),alpha('#065c69', 0.5)),pch=21,pt.cex=2.75,
       cex = 0.8,title = '', xpd = TRUE, inset = c(-0.21,0.1), bty = 'n')
# dev.off()
points(diet.aa$length, diet.aa$dbonds, cex = diet.aa$av.ko1*0.014+0.6, col = '#663333', bg = alpha('#663333', 0.5), pch = 21)
points(diet.aa$length, diet.aa$dbonds, cex = diet.aa$av.ko2*0.014+0.6, col = '#065c69', bg = alpha('#065c69', 0.5), pch = 21)


dev.off()






### INDIVIDUAL PLOT - PC_ae class
pdf('Distribution_of_PC_ae_class_DIET_AAV_INDIVID.pdf', w = 14, h = 5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0), mfrow=c(1,1)) 
par(mar = c(4,4,4,7), mgp = c(3,1,0), mfrow=c(2,3)) 

range(diet.ae$av.wt)
plot(diet.ae$length, diet.ae$dbonds, cex = diet.ae$av.wt*0.5+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'WT', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('A')), bty = "n", inset = c(-0.2, -0.2), xpd = TRUE)

axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','1','  5'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.05,1,5)*0.5+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
# dev.off()


# pdf('Distribution_of_PC_ae_class_DIET_AAV_INDIVID_KO.pdf', w = 6, h = 3.5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0),mfrow=c(1,1)) 
plot(diet.ae$length, diet.ae$dbonds, cex = diet.ae$av.ko1*0.5+0.6, col = '#663333', bg = alpha('#663333', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'KOS(1)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('B')), bty = "n", inset = c(-0.2, -0.2), xpd = TRUE)
axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','1','  5'),col='#663333',pt.bg = alpha('#663333', 0.5),
       pch=21,pt.cex=c(0.05,1,5)*0.5+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
# dev.off()


# pdf('Distribution_of_PC_ae_class_DIET_AAV_INDIVID_KOS.pdf', w = 6, h = 3.5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0),mfrow=c(1,1)) 
plot(diet.ae$length, diet.ae$dbonds, cex = diet.ae$av.ko2*0.5+0.6, col = '#065c69', bg = alpha('#065c69', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'KOS(2)', bty = "n", inset = c(-0.02, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('C')), bty = "n", inset = c(-0.2, -0.2), xpd = TRUE)
axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','1','  5'),col='#065c69',pt.bg = alpha('#065c69', 0.5),
       pch=21,pt.cex=c(0.05,1,5)*0.5+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
# dev.off()


par(mar = c(6,4,4,1))
b = barplot(diet.ae$av.wt/sum(diet.ae$av.wt)*100, col = alpha('#b4b4b4',0.8), yaxt = 'n', ylab = 'Percentage')
legend('topleft', legend = expression(bold('D')), bty = "n", inset = c(-0.15, -0.25), xpd = TRUE)
axis(2, las = 2)
axis(1, b, labels = rownames(diet.ae), las = 2)
legend('topleft', legend = 'WT', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)

b = barplot(diet.ae$av.ko1/sum(diet.ae$av.ko1)*100, col = alpha('#663333', 0.8), yaxt = 'n', ylab = 'Percentage')
legend('topleft', legend = expression(bold('E')), bty = "n", inset = c(-0.15, -0.25), xpd = TRUE)
axis(2, las = 2)
axis(1, b, labels = rownames(diet.ae), las = 2)
legend('topleft', legend = 'KOS(1)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)

b = barplot(diet.ae$av.ko2/sum(diet.ae$av.ko2)*100, col = alpha('#065c69', 0.8), yaxt = 'n', ylab = 'Percentage')
legend('topleft', legend = expression(bold('F')), bty = "n", inset = c(-0.15, -0.25), xpd = TRUE)
axis(2, las = 2)
axis(1, b, labels = rownames(diet.ae), las = 2)
legend('topleft', legend = 'KOS(2)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)

dev.off()


### OVERLAY PLOT - PC_ae class
pdf('Distribution_of_PC_ae_class_DIET_AAV_OVERLAY.pdf', w = 12, h = 2.5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0), mfrow=c(1,1)) 
par(mar = c(4,4,4,7), mgp = c(3,1,0), mfrow=c(1,3)) 

plot(diet.ae$length, diet.ae$dbonds, cex = diet.ae$av.wt*0.7+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'WT + KOS(1)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('A')), bty = "n", inset = c(-0.25, -0.2), xpd = TRUE)

axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','1','  5'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.05,1,5)*0.5+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
legend('bottomright', legend = paste0(c('WT', 'KOS(1)')),col=c('#b4b4b4', '#663333'),
       pt.bg = c(alpha('#b4b4b4', 0.5), alpha('#663333', 0.5)),pch=21,pt.cex=2.75,
       cex = 0.8,title = '', xpd = TRUE, inset = c(-0.19,0), bty = 'n')
# dev.off()
points(diet.ae$length, diet.ae$dbonds, cex = diet.ae$av.ko1*0.7+0.6, col = '#663333', bg = alpha('#663333', 0.5), pch = 21)

plot(diet.ae$length, diet.ae$dbonds, cex = diet.ae$av.wt*0.5+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'WT + KOS', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('B')), bty = "n", inset = c(-0.25, -0.2), xpd = TRUE)

axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','1','  5'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.05,1,5)*0.5+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
legend('bottomright', legend = paste0(c('WT', 'KOS(2)')),col=c('#b4b4b4', '#065c69'),
       pt.bg = c(alpha('#b4b4b4', 0.5), alpha('#065c69', 0.5)),pch=21,pt.cex=2.75,
       cex = 0.8,title = '', xpd = TRUE, inset = c(-0.21,0), bty = 'n')
# dev.off()
points(diet.ae$length, diet.ae$dbonds, cex = diet.ae$av.ko2*0.014+0.6, col = '#065c69', bg = alpha('#065c69', 0.5), pch = 21)


plot(diet.ae$length, diet.ae$dbonds, cex = diet.ae$av.wt*0.5+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,6.5))
legend('topleft', legend = 'WT + KOS(1) + KOS(2)', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('C')), bty = "n", inset = c(-0.25, -0.2), xpd = TRUE)

axis(1, seq(24,42,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.05','1','  5'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.05,1,5)*0.5+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
legend('bottomright', legend = paste0(c('WT', 'KOS(1)', 'KOS(2)')),col=c('#b4b4b4', '#663333','#065c69'),
       pt.bg = c(alpha('#b4b4b4', 0.5), alpha('#663333', 0.5),alpha('#065c69', 0.5)),pch=21,pt.cex=2.75,
       cex = 0.8,title = '', xpd = TRUE, inset = c(-0.21,0.1), bty = 'n')
# dev.off()
points(diet.ae$length, diet.ae$dbonds, cex = diet.ae$av.ko1*0.5+0.6, col = '#663333', bg = alpha('#663333', 0.5), pch = 21)
points(diet.ae$length, diet.ae$dbonds, cex = diet.ae$av.ko2*0.5+0.6, col = '#065c69', bg = alpha('#065c69', 0.5), pch = 21)

dev.off()