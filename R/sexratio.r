sexratio <- 
  function (species = "name species", cl = 1, unit = "cm", age = T, 
    n = 1000) 
{
    call <- match.call()
    cv <- function(data) {
        data <- as.numeric(data)
        md <- median(data)
        var <- (sum((data - md)^2))/(length(data) - 1)
        cv <- abs(var^0.5)/md
        cv[is.na(cv)] <- 0
        cv[is.infinite(cv)] <- NA
        res <- c(round(md, 3), round(cv, 3))
    }
    mcinf <- function(x) {
        k = 1
        x2 <- c()
        for (k in seq(along = x)) {
            if ((x[k] - floor(x[k])) < 0.5) 
                x[k] <- floor(x[k])
            else x[k] <- floor(x[k]) + 0.5
            x2 <- c(x2, x[k])
            k = k + 1
        }
        x2
    }
    freq.ori <- function(data) {
        t <- nrow(data[(data$sex == 1 | data$sex == 2), ])
        h <- nrow(data[data$sex == 2, ])
        freq <- round((h/t) * 100, 1)
        freq
    }
    
    freq.fun <- function(data, i) {
        dat <- data[i, ]
        dat <- as.data.frame(dat)
        freq.b <- freq.ori(dat)
        freq.b
    }
    dat <- read.csv(file = paste(species, 
        ".csv", sep = ""), header = T)
    dat <- data.frame(dat)
    sex.tal.dat <- data.frame(cbind(dat$tal, dat$sex))
    sex.tal.dat <- sex.tal.dat[sex.tal.dat[, 2] == 1 | sex.tal.dat[, 
        2] == 2, ]
    sex.tal.dat <- na.omit(sex.tal.dat)
    sex.tal.dat <- sex.tal.dat
    sex.tal.dat <- as.data.frame(sex.tal.dat)
    dimnames(sex.tal.dat)[[2]] <- c("tal", "sex")
    if (cl == 0.5) 
        sex.tal.dat$tal <- mcinf(sex.tal.dat$tal)
    if (cl == 1) 
        sex.tal.dat$tal <- floor(sex.tal.dat$tal)
    s <- unique(sex.tal.dat$tal)
    res <- c()
    rest <- matrix(ncol = 5)
    for (i in s) {
        d <- sex.tal.dat[sex.tal.dat[, 1] == i, ]
        d <- as.data.frame(d)
        N <- nrow(d)
        res.ori <- freq.ori(d)
        sex.boot <- boot(d, freq.fun, R = n)
        res.boot <- cv(sex.boot$t)
        i <- round(i, 1)
        res <- cbind(i, round(res.ori, 1), round(res.boot[1], 
            1), res.boot[2], N)
        rest <- matrix(rbind(rest, res), ncol = 5)
    }
    rest <- as.data.frame(rest)
    dimnames(rest)[[2]] <- c(paste("Length(", unid, ")",sep=""), "%original", 
        "%boot", "CV-boot", "n")
    rest <- rest[2:nrow(rest), ]
    ordt <- sort.list(rest[, 1])
    rest <- rest[ordt, ]
    rest.tot.perc.tal <- round((dim (subset(sex.tal.dat, sex.tal.dat$sex==2))[1]/ dim(sex.tal.dat)[1])*100,1)
    res.tot.tal <- as.data.frame(cbind(rest.tot.perc.tal, round(sum(rest$n * rest$"CV-boot")/sum(rest$n),4) , sum(rest$n)))
    colnames(res.tot.tal) <- c("%original","CV(weihgted mean)", "n")
    
    if (edad == T) {
        sex.eda.dat <- data.frame(cbind(dat$eda, dat$sex))
        sex.eda.dat <- sex.eda.dat[sex.eda.dat[, 2] == 1 | sex.eda.dat[, 
            2] == 2, ]
        sex.eda.dat <- na.omit(sex.eda.dat)
        sex.eda.dat <- floor(sex.eda.dat)
        sex.eda.dat <- as.data.frame(sex.eda.dat)
        dimnames(sex.eda.dat)[[2]] <- c("eda", "sex")
        e <- unique(sex.eda.dat$eda)
        res <- c()
        rese <- matrix(ncol = 5)
        for (i in e) {
            d <- sex.eda.dat[sex.eda.dat[, 1] == i, ]
            d <- as.data.frame(d)
            N <- nrow(d)
            rese.ori <- freq.ori(d)
            sex.e.boot <- boot(d, freq.fun, R = n)
            rese.boot <- cv(sex.e.boot$t)
            i <- round(i, 0)
            res <- cbind(i, round(rese.ori, 1), round(rese.boot[1], 
                1), rese.boot[2], N)
            rese <- matrix(rbind(rese, res), ncol = 5)
        }
        rese <- as.data.frame(rese)
        dimnames(rese)[[2]] <- c("age", "%original", " %boot", 
            "CV-boot", "n")
        rese <- rese[2:nrow(rese), ]
        orde <- sort.list(rese[, 1])
        rese <- rese[orde, ]
    }
    res.tot.perc.eda <- round((dim (subset(sex.eda.dat, sex.eda.dat$sex==2))[1]/ dim(sex.eda.dat)[1])*100,1)
    res.tot.eda <- as.data.frame(cbind(res.tot.perc.eda, round(sum(rese$n * rese$"CV-boot", na.rm=T)/sum(rese$n),4), sum(rese$n)))
    colnames(res.tot.eda) <- c("%original","CV(weighted mean)", "n")
    path <- getwd()
    sink(file = paste(species, "_sexratio_results.txt", sep = ""))
    cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", "Time:", 
        format(Sys.time(), "%X"))
    cat("\n")
    cat("\n")
    cat("\n")
    cat("SEXRATIO (%Females)", "\n")
    cat("\n")
    cat("**********************************************************", 
        "\n")
    cat("species:", species, "\n")
    cat("\n")
    cat("Arguments Routine:", "\n")
    print(call)
    cat("**********************************************************", 
        "\n")
    cat("\n")
    cat("\n", "SEXRATIO BY SIZE", "\n")
    cat("\n")
    print(rest, row.names=FALSE)
    cat("\n")
    cat("\n", "GLOBAL:", "\n")
    print(res.tot.tal, row.names=FALSE)
    cat("\n")
    cat("\n")
    
    if (edad == T) {
        cat("\n", "SEXRATIO BY AGE", "\n")
        cat("\n")
        print(rese, row.names=FALSE)
        cat("\n")
        cat("\n", "GLOBAL:", "\n")
        print(res.tot.eda, row.names=FALSE)
    }
    sink()
    pdf(file = paste(species,"_sexratio_plots.pdf", sep = ""), width = 9, height = 13)
    layout(matrix(c(1, 2), 2, 1, byrow = T), widths = c(3, 3), heights = c(1.5, 1.5))    
    par(mar=c(5,4,4,5)+.1, oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
        cex.lab = 1.2)
    rest.plot <- ifelse (rest[,4]==0, 0.001, rest[,4])
    plot(rest.plot ~ rest[, 1], type="h",col="red",xaxt="n",yaxt="n",xlab="",ylab="")
    axis(4)
    mtext("CV",side=4,line=3)
    par(new=TRUE)
    plot(rest[, 2] ~ rest[, 1], xlab = paste("Talla (",unid, ")",sep=""), ylab = " % Hembras", main = "Sexratio Talla",type="b", lty=1,pch=19, col="navy")
    legend("top",col=c("navy","red"),lty=1,pch=c(19,-1), merge=TRUE, bty="n", legend=c("% Hembras","CV"))
    if (edad == T) {
      rese.plot <- ifelse (rese[,4]==0, 0.002, rese[,4])
      plot(rese.plot ~ rese[, 1], type="h",col="red",ylim=c(0, max(rese[,4])*1.08),xaxt="n",yaxt="n",xlab="",ylab="")
      axis(4)
      mtext("CV",side=4,line=3)
      par(new=TRUE)
      plot(rese[, 2] ~ rese[, 1], xlab = "Edad (años)", ylab = " % Hembras", 
            main = "Sexratio Edad", type="b", pch=19, col="darkgreen")
      legend("top",col=c("darkgreen","red"),lty=1,pch=c(19,-1), merge=TRUE, bty="n",legend=c("% Hembras","CV"))
            }
  
    graphics.off()
    cat("\n")
    cat(paste("EL PROGRAMA HA TERMINADO.\n","PUEDES VER LOS RESULTADOS EN EL DIRECTORIO DE TRABAJO:", path))
    cat("\n")
}

