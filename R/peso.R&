"peso" <-
function (especie = "nombre especie", cl = 1, unid = "cm", sex = F, 
    a = 5e-04, b = 3, n = 1000) 
{
    call <- match.call()
    cv <- function(data) {
        data <- as.numeric(data)
        med <- round(median(data), 6)
        var <- (sum((data - med)^2))/(length(data) - 1)
        cv <- round(abs(var^0.5/med), 3)
        list(mediana = med, cv = cv)
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
    path <- getwd()
    plot.ib <- function(x, which = 1:2, main = main, panel = points, 
        ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
        ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75) {
        show <- rep(FALSE, 2)
        show[which] <- TRUE
        r <- residuals(x)
        r <- (r - mean(r))/sd(r)
        yh <- predict(x)
        n <- length(r)
        if (id.n > 0) {
            if (is.null(labels.id)) 
                labels.id <- paste(1:n)
            iid <- 1:id.n
            show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
            if (any(show[2])) 
                show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
            text.id <- function(x, y, ind, adj.x = FALSE) text(x - 
                if (adj.x) 
                  strwidth(" ") * cex.id
                else 0, y, labels.id[ind], cex = cex.id, xpd = TRUE, 
                adj = if (adj.x) 
                  1)
        }
        one.fig <- prod(par("mfcol")) == 1
        if (ask) {
            op <- par(ask = TRUE)
            on.exit(par(op))
        }
        if (show[1]) {
            ylim <- range(r, na.rm = TRUE)
            if (id.n > 0) 
                ylim <- ylim + c(-1, 1) * 0.08 * diff(ylim)
            plot(yh, r, xlab = "Predicciones", ylab = "Resíduos estandarizados", 
                main = main, ylim = ylim, type = "n", ...)
            panel(yh, r, ...)
            if (id.n > 0) {
                y.id <- r[show.r]
                y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
                text.id(yh[show.r], y.id, show.r, adj.x = TRUE)
            }
            abline(h = 0, lty = 3, col = "red")
        }
        if (show[2]) {
            ylim <- range(r, na.rm = TRUE)
            ylim[2] <- ylim[2] + diff(ylim) * 0.075
            qq <- qqnorm(r, ylab = "Resíduos estandarizados", 
                xlab = "Cuantiles teóricos", main = main, ylim = ylim, 
                ...)
            qqline(r, lty = 3)
            if (id.n > 0) 
                text.id(qq$x[show.r], qq$y[show.r], show.r, adj.x = TRUE)
        }
        invisible()
    }
    pes.f <- function(data) {
        pes.dat <- data.frame(cbind(data$tal, data$pes))
        pes.dat <- na.omit(pes.dat)
        pes.dat <- as.data.frame(pes.dat)
        n <- nrow(pes.dat)
        dimnames(pes.dat)[[2]] <- c("tal", "pes")
        attach(pes.dat)
        pes.ori <- nls(pes ~ a * tal^b, start = list(a = a, b = b))
        pes.fit <- fitted(pes.ori)
        pes.fit <- as.vector(pes.fit)
        pes.fun <- function(data) coef(nls(pes ~ a * tal^b, start = coef(pes.ori), 
            data = data))
        pes.case <- function(data, i) pes.fun(data[i, ])
        pes.boot <- boot(pes.dat, pes.case, R = n)
        a <- pes.boot$t[, 1]
        b <- pes.boot$t[, 2]
        a.boot <- cv(a)
        b.boot <- cv(b)
        par.boot <- rbind(a.boot, b.boot)
        par.ori <- round(coef(pes.ori), 6)
        par <- cbind(par.ori, par.boot)
        dimnames(par)[[1]] <- c("ordenada", "pendiente")
        dimnames(par)[[2]] <- c("Estima datos", "Estima boot", 
            "CV boot")
        detach(pes.dat)
        res <- list(n = n, dat = pes.dat, nls.ori = pes.ori, 
            nls.boot = pes.boot, parametros = par)
        rm(pes.dat, pes.ori, pes.fit, pes.boot, pes.case, a, 
            b, a.boot, b.boot, par.ori, par.boot)
        res
    }
    dat <- read.csv(file = paste(especie,".csv", sep = ""), header = T)
    dat <- data.frame(dat)
    if (cl == 0.5) 
        dat$tal <- mcinf(dat$tal)
    if (cl == 1) 
        dat$tal <- floor(dat$tal)
    dat$tal <- dat$tal + (cl/2)
    orde <- sort.list(dat[, 1])
    dat <- dat[orde, ]
    if (sex == F) {
        pes.res <- pes.f(dat)
        pdf(file = paste(especie, "_peso_plots.pdf", sep = ""), width = 9, height = 13)
        layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 4, 2, byrow = T), 
            widths = c(3, 3), heights = c(1.5, 1.5))
        par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
            cex.lab = 1.2)
        plot(pes.res$dat$pes ~ pes.res$dat$tal, xlab = paste("Talla (", 
            unid, ")",sep=""), ylab = " Peso (g)", main = "Relación Talla-Peso", 
            pch = "+", col="darkgrey", ylim=c(min(pes.res$dat$pes)*0.9, max(pes.res$dat$pes)*1.1))
        lines(pes.res$dat$tal, fitted(pes.res$nls.ori), type = "l", 
            col = "red")
        res.acf <- residuals(pes.res$nls.ori)
        x.acf <- sample(res.acf, length(res.acf), replace = F)
        acf(x.acf, main = "Correlograma")
        plot.ib(pes.res$nls.ori, which = 1, main = "Resíduos vs Predicciones")
        plot.ib(pes.res$nls.ori, which = 2, main = "Normal QQ Plot")
        hist(pes.res$nls.boot$t[, 1], main = "BOOTSTRAP", xlab = "a", 
            ylab = "Frecuencia", col = "grey")
        hist(pes.res$nls.boot$t[, 2], main = "BOOTSTRAP", xlab = "b", 
            ylab = "Frecuencia", col = "grey")
        qqnorm(pes.res$nls.boot$t[, 1], xlab = "Cuantiles teóricos", 
            ylab = "a", main = "")
        qqline(pes.res$nls.boot$t[, 1])
        qqnorm(pes.res$nls.boot$t[, 2], xlab = "Cuantiles teóricos", 
            ylab = "b", main = "")
        qqline(pes.res$nls.boot$t[, 2])
        mtext(outer = T, "Relación Talla-Peso\nGráficos Diagnóstico" , 
            side = 3, cex = 1.2)
        dev.off()
        sink(file = paste(especie, "_peso_resultados.txt", sep = ""))
        cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", "Hora:", 
            format(Sys.time(), "%X"))
        cat("\n")
        cat("\n")
        cat("\n")
        cat("RELACIÓN TALLA-PESO: AJUSTE NO-LINEAL", "\n")
        cat("\n")
        cat("*******************************************************", 
            "\n")
        cat("Especie:", especie, "\n")
        cat("Sexo: Ambos", "\n")
        cat("N:", pes.res$n, "\n")
        cat("\n")
        cat("Argumentos Rutina:", "\n")
        print(call)
        cat("*******************************************************", 
            "\n")
        cat("\n")
        cat("RESULTADOS AJUSTE (datos originales) ", "\n")
        cat("\n")
        print(pes.res$nls.ori)
        cat("\n")
        cat("\n")
        cat("\n")
        cat("\n")
        cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP", "\n")
        cat("\n")
        cat("\n")
        print(pes.res$parametros)
        sink()
        rm(pes.res)
    }
    if (sex == T) {
        datm <- dat[dat[, 5] == 1, ]
        datf <- dat[dat[, 5] == 2, ]
        rm(dat)
        pes.m.res <- pes.f(datm)
        pdf(file = paste(especie, "_pesomachos_plots.pdf", sep = ""), width = 9, 
            height = 13)
        layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 4, 2, byrow = T), 
            widths = c(3, 3), heights = c(1.5, 1.5))
        par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
            cex.lab = 1.2)
        plot(pes.m.res$dat$pes ~ pes.m.res$dat$tal, xlab = paste("Talla (", 
            unid, ")", sep=""), ylab = " Peso (g)", main = "Relación Talla-Peso", 
            pch = "+", col="darkgrey",ylim=c(min(pes.m.res$dat$pes)*0.9, max(pes.m.res$dat$pes)*1.1))
        lines(pes.m.res$dat$tal, fitted(pes.m.res$nls.ori), type = "l", 
            col = "red")
        res.acf <- residuals(pes.m.res$nls.ori)
        x.acf <- sample(res.acf, length(res.acf), replace = F)
        acf(x.acf, main = "Correlograma")
        plot.ib(pes.m.res$nls.ori, which = 1, main = "Resíduos vs Predicciones")
        plot.ib(pes.m.res$nls.ori, which = 2, main = "Normal QQ Plot")
        hist(pes.m.res$nls.boot$t[, 1], main = "BOOTSTRAP", xlab = "a", 
            ylab = "Frecuencia", col = "blue")
        hist(pes.m.res$nls.boot$t[, 2], main = "BOOTSTRAP", xlab = "b", 
            ylab = "Frecuencia", col = "blue")
        qqnorm(pes.m.res$nls.boot$t[, 1], xlab = "Cuantiles teóricos", 
            ylab = "a", main = "")
        qqline(pes.m.res$nls.boot$t[, 1])
        qqnorm(pes.m.res$nls.boot$t[, 2], xlab = "Cuantiles teóricos", 
            ylab = "b", main = "")
        qqline(pes.m.res$nls.boot$t[, 2])
        mtext(outer = T, "MACHOS: Relación Talla-Peso\nGráficos Diagnóstico" , 
            side = 3, cex = 1.2)
        dev.off()
        sink(file = paste(especie, "_pesomachos_resultados.txt", sep = ""))
        cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", "Hora:", 
            format(Sys.time(), "%X"))
        cat("\n")
        cat("\n")
        cat("\n")
        cat("RELACIÓN TALLA-PESO: AJUSTE NO-LINEAL", "\n")
        cat("\n")
        cat("*******************************************************", 
            "\n")
        cat("Especie:", especie, "\n")
        cat("Sexo: Machos", "\n")
        cat("N:", pes.m.res$n, "\n")
        cat("\n")
        cat("Argumentos Rutina:", "\n")
        print(call)
        cat("*******************************************************", 
            "\n")
        cat("\n")
        cat("RESULTADOS AJUSTE (datos originales) ", "\n")
        cat("\n")
        print(pes.m.res$nls.ori)
        cat("\n")
        cat("\n")
        cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP ", "\n")
        cat("\n")
        print(pes.m.res$parametros)
        sink()
        rm(pes.m.res)
        pes.f.res <- pes.f(datf)
        pdf(file = paste(especie, "_pesohembras_plots.pdf", sep = ""), width = 9, 
            height = 13)
        layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 4, 2, byrow = T), 
            widths = c(3, 3), heights = c(1.5, 1.5))
        par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
            cex.lab = 1.2)
        plot(pes.f.res$dat$pes ~ pes.f.res$dat$tal, xlab = paste("Talla (", 
            unid, ")",sep=""), ylab = " Peso (g)", main = "Relación Talla-Peso", 
            pch = "+",col="darkgrey",ylim=c(min(pes.f.res$dat$pes)*0.9, max(pes.f.res$dat$pes)*1.1))
        lines(pes.f.res$dat$tal, fitted(pes.f.res$nls.ori), type = "l", 
            col = "red")
        res.acf <- residuals(pes.f.res$nls.ori)
        x.acf <- sample(res.acf, length(res.acf), replace = F)
        acf(x.acf, main = "Correlograma")
        plot.ib(pes.f.res$nls.ori, which = 1, main = "Resíduos vs Predicciones")
        plot.ib(pes.f.res$nls.ori, which = 2, main = "Normal QQ Plot")
        hist(pes.f.res$nls.boot$t[, 1], main = "BOOTSTRAP", xlab = "a", 
            ylab = "Frecuencia", col = "pink")
        hist(pes.f.res$nls.boot$t[, 2], main = "BOOTSTRAP", xlab = "b", 
            ylab = "Frecuencia", col = "pink")
        qqnorm(pes.f.res$nls.boot$t[, 1], xlab = "Cuantiles teóricos", 
            ylab = "a", main = "")
        qqline(pes.f.res$nls.boot$t[, 1])
        qqnorm(pes.f.res$nls.boot$t[, 2], xlab = "Cuantiles teóricos", 
            ylab = "b", main = "")
        qqline(pes.f.res$nls.boot$t[, 2])
        mtext(outer = T, "HEMBRAS: Relación Talla-Peso\nGráficos Diagnóstico" , 
            side = 3, cex = 1.2)
        dev.off()
        sink(file = paste(especie, "_pesohembras_resultados.txt", sep = ""))
        cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", "Hora:", 
            format(Sys.time(), "%X"))
        cat("\n")
        cat("\n")
        cat("\n")
        cat("RELACIÓN TALLA-PESO: AJUSTE NO-LINEAL", "\n")
        cat("\n")
        cat("*******************************************************", 
            "\n")
        cat("Especie:", especie, "\n")
        cat("Sexo: Hembras", "\n")
        cat("N:", pes.f.res$n, "\n")
        cat("\n")
        cat("Argumentos Rutina:", "\n")
        print(call)
        cat("*******************************************************", 
            "\n")
        cat("\n")
        cat("RESULTADOS AJUSTE (datos originales) ", "\n")
        cat("\n")
        print(pes.f.res$nls.ori)
        cat("\n")
        cat("\n")
        cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP ", "\n")
        cat("\n")
        print(pes.f.res$parametros)
        sink()
        rm(pes.f.res)
    }
    graphics.off()
    cat("\n")
    cat(paste("EL PROGRAMA HA TERMINADO.\n","PUEDES VER LOS RESULTADOS EN EL DIRECTORIO DE TRABAJO:", path))
    cat("\n")
}

