#examine the difference in raw counts between tumor and non-tumor
#row 6648 in merged.table = ENSG00000132698 = rab25
rab25 <- merged.tables[6648,]
dev.new()
par(mar=c(12, 5, 4, 2))
matplot(t(rab25), type="p", pch=16, col=1, main="rab25 raw counts", ylab="raw counts", xaxt="n", log="y")
axis(1, at=1:24, labels=colnames(rab25), las=2)