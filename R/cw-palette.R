tol5qualitative= c("#332288","#88CCEE","#117733","#DDCC77","#CC6677")
tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")

redblue<-colorRampPalette(c("red","white","blue"))
    magblue<-colorRampPalette(c("magenta","white","blue"))
    ## grnvi<-colorRampPalette(c("#00441b","gray90","#40004b"))
grnmag<-colorRampPalette(c("#00441b","white","magenta"))

brn <- "#543005"; grn <- "#003c30"
brngrn <- colorRampPalette(c("#543005","white","#003c30"))

lighten <- function(color, factor=1.2){
    col <- col2rgb(color)
    col <- col*factor
    col <- rgb(t(col), maxColorValue=255)
    col
}


darken <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}

myblue <- "#5750a2" # from olly
myred <- "#842525" # from olly
## mygreen <- "#2e8539" ##00441b"
mygreen= darken("#8fd175")
grnvi<-colorRampPalette(c(darken(myblue),"gray95",darken(myred))) # updated to match Olly's figure 1
grnvi<-colorRampPalette(c(lighten(myblue),"gray95",lighten(myred))) # updated to match Olly's figure 1
## grnvi<-colorRampPalette(c(myblue,"gray95",myred)) # updated to match Olly's figure 1
