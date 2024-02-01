library(readr)
library(scales)

diet = read.delim('../diet_car_concentrations.csv', stringsAsFactors = F, sep = ',')
ko.colnames = colnames(diet)[grep('KO', colnames(diet))]
kos.colnames = colnames(diet)[grep('KOS', colnames(diet))]
ko.colnames = setdiff(ko.colnames, kos.colnames)

for (i in 1:nrow(diet)){
  diet$av.wt[i] =  mean(as.numeric(diet[i,grep('WT',colnames(diet))]))
  diet$av.ko1[i] =  mean(as.numeric(diet[i,ko.colnames]))
  diet$av.ko2[i] =  mean(as.numeric(diet[i,kos.colnames]))
}

diet$length = gsub('C', '', sapply(strsplit(rownames(diet), ' '), '[', 1) )                      
diet$dbonds = sapply(strsplit(diet$length, ':'), '[', 2)       
diet$dbonds[c(1,2,3,6,7,8,11,13)] = 0
diet$dbonds[10] = NA
diet$length = sapply(strsplit(diet$length, ':'), '[', 1)                      
diet$length[c(4,5)] = NA
diet$length[10] = NA

### INDIVIDUAL PLOT - CAR class
pdf('Distribution_of_CAR_class_DIET_INDIVID.pdf', w = 14, h = 5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0), mfrow=c(1,1)) 
par(mar = c(4,4,4,7), mgp = c(3,1,0), mfrow=c(2,3)) 

plot(diet$length, diet$dbonds, cex = diet$av.wt*0.4+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length',
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,2.5), xlim = c(-1, 18.5))
legend('topleft', legend = 'WT', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('A')), bty = "n", inset = c(-0.2, -0.2), xpd = TRUE)

axis(1, seq(0,40,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.1','1',' 5'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.1,1,5)*0.4+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
# dev.off()
range(diet$av.wt)

# pdf('Distribution_of_CAR_class_DIET_INDIVID_KO.pdf', w = 6, h = 3.5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0), mfrow=c(1,1)) 
plot(diet$length, diet$dbonds, cex = diet$av.ko1*0.4+0.6, col = '#663333', bg = alpha('#663333', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,2.5), xlim = c(-1, 18.5))
legend('topleft', legend = 'KO', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('B')), bty = "n", inset = c(-0.2, -0.2), xpd = TRUE)
axis(1, seq(0,40,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.1','1',' 5'),col='#663333',pt.bg = alpha('#663333', 0.5),
       pch=21,pt.cex=c(0.1,1,5)*0.4+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
# dev.off()


# pdf('Distribution_of_CAR_class_DIET_INDIVID_KOS.pdf', w = 6, h = 3.5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0), mfrow=c(1,1)) 
plot(diet$length, diet$dbonds, cex = diet$av.ko2*0.4+0.6, col = '#065c69', bg = alpha('#065c69', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,2.5), xlim = c(-1, 18.5))
legend('topleft', legend = 'KOS', bty = "n", inset = c(-0.02, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('C')), bty = "n", inset = c(-0.2, -0.2), xpd = TRUE)
axis(1, seq(0,40,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.1','1',' 5'),col='#065c69',pt.bg = alpha('#065c69', 0.5),
       pch=21,pt.cex=c(0.01,1,5)*0.4+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
# dev.off()


par(mar = c(8,4,4,1))
b = barplot(diet$av.wt/sum(diet$av.wt)*100, col = alpha('#b4b4b4',0.8), yaxt = 'n', ylab = 'Percentage')#, ylim = c(0,30))
legend('topleft', legend = expression(bold('D')), bty = "n", inset = c(-0.15, -0.25), xpd = TRUE)
axis(2, las = 2)
axis(1, b, labels = rownames(diet), las = 2)
legend('topleft', legend = 'WT', bty = "n", inset = c(0.12, 0.0), xpd = TRUE)

b = barplot(diet$av.ko1/sum(diet$av.ko1)*100, col = alpha('#663333', 0.8), yaxt = 'n', ylab = 'Percentage')#, ylim = c(0,30))
legend('topleft', legend = expression(bold('E')), bty = "n", inset = c(-0.15, -0.25), xpd = TRUE)
axis(2, las = 2)
axis(1, b, labels = rownames(diet), las = 2)
legend('topleft', legend = 'KO', bty = "n", inset = c(0.12, 0.0), xpd = TRUE)

b = barplot(diet$av.ko2/sum(diet$av.ko2)*100, col = alpha('#065c69', 0.8), yaxt = 'n', ylab = 'Percentage')#, ylim = c(0,30))
legend('topleft', legend = expression(bold('F')), bty = "n", inset = c(-0.15, -0.25), xpd = TRUE)
axis(2, las = 2)
axis(1, b, labels = rownames(diet), las = 2)
legend('topleft', legend = 'KOS', bty = "n", inset = c(0.12, 0.0), xpd = TRUE)

dev.off()



### OVERLAY PLOT - CAR class
pdf('Distribution_of_CAR_class_DIET_OVERLAY.pdf', w = 12, h = 2.5)
# par(mar = c(4,4,4,6), mgp = c(3,1,0), mfrow=c(1,1)) 
par(mar = c(4,4,4,7), mgp = c(3,1,0), mfrow=c(1,3)) 

plot(diet$length, diet$dbonds, cex = diet$av.wt*0.4+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,2.5), xlim = c(-1, 18.5))
legend('topleft', legend = 'WT + KO', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('A')), bty = "n", inset = c(-0.25, -0.2), xpd = TRUE)

axis(1, seq(0,40,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.1','1',' 5'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.1,1,5)*0.4+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
legend('bottomright', legend = paste0(c('WT', 'KO')),col=c('#b4b4b4', '#663333'),
       pt.bg = c(alpha('#b4b4b4', 0.5), alpha('#663333', 0.5)),pch=21,pt.cex=2.75,
       cex = 0.8,title = '', xpd = TRUE, inset = c(-0.19,0), bty = 'n')
# dev.off()
points(diet$length, diet$dbonds, cex = diet$av.ko1*0.4+0.6, col = '#663333', bg = alpha('#663333', 0.5), pch = 21)

plot(diet$length, diet$dbonds, cex = diet$av.wt*0.4+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,2.5), xlim = c(-1, 18.5))
legend('topleft', legend = 'WT + KOS', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('B')), bty = "n", inset = c(-0.25, -0.2), xpd = TRUE)

axis(1, seq(0,40,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.1','1',' 5'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.1,1,5)*0.4+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
legend('bottomright', legend = paste0(c('WT', 'KOS')),col=c('#b4b4b4', '#065c69'),
       pt.bg = c(alpha('#b4b4b4', 0.5), alpha('#065c69', 0.5)),pch=21,pt.cex=2.75,
       cex = 0.8,title = '', xpd = TRUE, inset = c(-0.21,0), bty = 'n')
# dev.off()
points(diet$length, diet$dbonds, cex = diet$av.ko2*0.4+0.6, col = '#065c69', bg = alpha('#065c69', 0.5), pch = 21)


plot(diet$length, diet$dbonds, cex = diet$av.wt*0.4+0.6, col = '#b4b4b4', bg = alpha('#b4b4b4', 0.5), pch = 21,
     main = '', ylab = 'Double bonds', xlab = 'Carbon chain length', 
     yaxt = 'n', xaxt = 'n', ylim = c(-0.5,2.5), xlim = c(-1, 18.5))
legend('topleft', legend = 'WT + KO + KOS', bty = "n", inset = c(-0.01, 0.0), xpd = TRUE)
legend('topleft', legend = expression(bold('C')), bty = "n", inset = c(-0.25, -0.2), xpd = TRUE)

axis(1, seq(0,40,2))
axis(2,seq(0,6,1), las = 2)
legend('topright', legend = c('0.1','1',' 5'),col='#b4b4b4',pt.bg = alpha('#b4b4b4', 0.5),
       pch=21,pt.cex=c(0.1,1,5)*0.4+0.6,cex = 0.8,title = 'Concentration', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
legend('bottomright', legend = paste0(c('WT', 'KO', 'KOS')),col=c('#b4b4b4', '#663333','#065c69'),
       pt.bg = c(alpha('#b4b4b4', 0.5), alpha('#663333', 0.5),alpha('#065c69', 0.5)),pch=21,pt.cex=2.75,
       cex = 0.8,title = '', xpd = TRUE, inset = c(-0.21,0.1), bty = 'n')
# dev.off()
points(diet$length, diet$dbonds, cex = diet$av.ko1*0.4+0.6, col = '#663333', bg = alpha('#663333', 0.5), pch = 21)
points(diet$length, diet$dbonds, cex = diet$av.ko2*0.4+0.6, col = '#065c69', bg = alpha('#065c69', 0.5), pch = 21)

dev.off()