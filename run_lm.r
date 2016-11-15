source("./script/func.regression.r");
args<-commandArgs(T)
# print(args)

ex.file=args[2]  # expression.file
me.file=args[3]  # methylation.file
# ped.file=commandArgs()[5]   # pedigree.file
# out.file=commandArgs()[6]   # output.file

write(paste("ex.file=", ex.file), stdout())
write(paste("me.file=", me.file), stdout())

obj <- func.glmnet.1.cv(ex.file, me.file);
