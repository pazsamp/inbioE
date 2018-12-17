"crecimiento" <-
function (especie = "nombre especie", cl = 1, unid = "cm", sex = F, 
    b = 3, Li = 80, Lfija = F, Ki = 0.7, T0i = -1, Wi = 2500, 
    Wfijo = F, Kwi = 0.4, T0wi = -1, n = 1000) 
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
    cre.tal.f <- function(data, sexo) {
      cre.tal.dat <- data.frame(cbind(data$tal, data$eda))
      cre.tal.dat <- na.omit(cre.tal.dat)
      cre.tal.dat <- as.data.frame(cre.tal.dat)
      ntal <- nrow(cre.tal.dat)
      dimnames(cre.tal.dat)[[2]] <- c("tal", "eda")
      attach(cre.tal.dat)
      if (Lfija == F) {
        cre.tal.ori <- nls(tal ~ L * (1 - exp(-K * (eda - 
                                                      T0))), start = list(L = Li, K = Ki, T0 = T0i), 
                           control = nls.control(maxiter = 100, tol = 1e-04))
        detach(cre.tal.dat)
        cre.tal.fit <- coef(cre.tal.ori)
        cre.tal.fit <- as.vector(cre.tal.fit)
        names(cre.tal.fit) <- c(paste("Linf(", unid, ")",sep=""), 
                                "K(año-1)", "t0(año)")
        Li.ori <- cre.tal.fit[1]
        K.ori <- cre.tal.fit[2]
        t0.ori <- cre.tal.fit[3]
        cre.tal.fun <- function(data) coef(nls(tal ~ L * 
                                                 (1 - exp(-K * (eda - T0))), start = coef(cre.tal.ori), 
                                               control = nls.control(maxiter = 100, tol = 1e-04), 
                                               data = data))
        cre.tal.case <- function(data, i) cre.tal.fun(data[i,])
        cre.tal.boot <- boot(cre.tal.dat, cre.tal.case, R = n)
        Li.b <- cre.tal.boot$t[, 1]
        K.tal.b <- cre.tal.boot$t[, 2]
        t0.b <- cre.tal.boot$t[, 3]
        Li.boot <- as.vector(cv(Li.b))
        K.tal.boot <- as.vector(cv(K.tal.b))
        t0.boot <- as.vector(cv(t0.b))
        names(Li.boot) <- c("Li: estima bootstrap", "Li:cv")
        names(K.tal.boot) <- c("K: estima bootstrap", "K:cv")
        names(t0.boot) <- c("t0: estima bootstrap", "t0:cv")
        par.boot <- rbind(Li.boot, K.tal.boot, t0.boot)
        par.ori <- round(coef(cre.tal.ori), 6)
        par.tal <- cbind(par.ori, par.boot)
        dimnames(par.tal)[[1]] <- c(paste("Linf(", unid, ")",sep=""), "k (año-1)", "t0 (año)")
        dimnames(par.tal)[[2]] <- c("Estima original", "Estima boot", "CV boot")
        if (sexo == "Ambos") 
          pdf(file = paste(especie,"_cretal_plots.pdf", sep = ""), width = 9, 
                  height = 13)
        if (sexo == "Machos") 
            pdf(file = paste(especie, "_cretalmachos_plots.pdf", sep = ""), width = 9, 
                  height = 13)
        if (sexo == "Hembras") 
            pdf(file = paste(especie, "_cretalhembras_plots.pdf", sep = ""), width = 9, height = 13)
          
        layout(matrix(c(1, 1, 2, 3, 4, 0), 3, 2, byrow = T), 
                widths = c(3, 3), heights = c(1.5, 1.5))
        par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
                cex.lab = 1.2)
        plot(cre.tal.dat$tal ~ cre.tal.dat$eda, xlab = "Edad (años)", ylab = paste("Talla (", unid, ")",sep=""), main = "Crecimiento en Talla", 
                pch = 1, ylim=c(min(cre.tal.dat$tal)*0.9,max(cre.tal.dat$tal)*1.1))
        lines(cre.tal.dat$eda, cre.tal.fit[1] * (1 - exp(-cre.tal.fit[2] * (cre.tal.dat$eda - cre.tal.fit[3]))), type = "l", 
            col = "red", lty = 2, lwd = 2)
        lines(cre.tal.dat$eda, median(Li.b) * (1 - exp(-median(K.tal.b) * (cre.tal.dat$eda - median(t0.b)))), type = "l", 
                col = "blue", lty = 3, lwd = 2)
        legend("bottomright",c("Predicción Original", "Predicción Bootstrap"), lwd = 1, lty = 1, col = c("red", "blue"), cex=0.8, bty="n")
            mtext(outer = T, " Crecimiento en Talla: Gráficos Diagnóstico", side = 3, cex = 1.2)     
        
        res.acf <- residuals(cre.tal.ori)
        x.acf <- sample(res.acf, length(res.acf), replace = F)
        acf(x.acf, main = "Correlograma")
        plot.ib(cre.tal.ori, which = 1, main = "Resíduos vs Predicciones")
        plot.ib(cre.tal.ori, which = 2, main = "Normal QQ Plot")
        layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = F),widths = c(3, 3), heights = c(1.5, 1.5))
        par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, cex.lab = 1.2)
        hist(Li.b, main = "BOOTSTRAP", xlab = paste("Linf (", unid, ")",sep=""), ylab = "Frecuencia")
        hist(K.tal.b, main = "BOOTSTRAP", xlab = "K (año-1)", ylab = "Frecuencia")
        hist(t0.b, main = "BOOTSTRAP", xlab = "t0 (año)", ylab = "Frecuencia")
        qqnorm(Li.b, xlab = "Cuantiles teóricos", ylab = "Cuantiles-Linf ", main = "BOOTSTRAP")
        qqline(Li.b)
        qqnorm(K.tal.b, xlab = "Cuantiles teóricos", ylab = "Cuantiles-K ", 
                main = "BOOTSTRAP")
        qqline(K.tal.b)
        qqnorm(t0.b, xlab = "Cuantiles teóricos", ylab = "Cuantiles-t0 ", 
                main = "BOOTSTRAP")
        qqline(t0.b)
        mtext(outer = T, " Crecimiento en Talla: Gráficos Diagnóstico" , 
                side = 3, cex = 1.2)
        dev.off()
        }
      if (Lfija == T) {
        attach(cre.tal.dat)
        cre.tal.ori <- nls(tal ~ Li * (1 - exp(-K * (eda - T0))), start = list(K = Ki, T0 = T0i), control = nls.control(maxiter = 100, tol = 1e-04))
        detach(cre.tal.dat)
        cre.tal.fit <- coef(cre.tal.ori)
        cre.tal.fit <- as.vector(cre.tal.fit)
        names(cre.tal.fit) <- c("K (año-1)", "t0 (año)")
        K.ori <- cre.tal.fit[1]
        t0.ori <- cre.tal.fit[2]
        cre.tal.fun <- function(data) coef(nls(tal ~ 
                                                 Li * (1 - exp(-K * (eda - T0))), start = coef(cre.tal.ori), 
                                               control = nls.control(maxiter = 100, tol = 1e-04), data=data))
        cre.tal.case <- function(data, i) cre.tal.fun(data[i,])
        cre.tal.boot <- boot(cre.tal.dat, cre.tal.case, R = n)
        
        K.tal.b <- cre.tal.boot$t[, 1]
        t0.b <- cre.tal.boot$t[, 2]
        K.tal.boot <- as.vector(cv(K.tal.b))
        t0.boot <- as.vector(cv(t0.b))
        names(K.tal.boot) <- c("K: estima bootstrap", "K:cv")
        names(t0.boot) <- c("t0: estima bootstrap", "t0:cv")
        par.boot <- rbind(K.tal.boot, t0.boot)
        par.ori <- round(coef(cre.tal.ori), 6)
        par.tal <- cbind(par.ori, par.boot)
        dimnames(par.tal)[[1]] <- c("k(año-1)", "t0(año)")
        dimnames(par.tal)[[2]] <- c("Estima original", "Estima boot", "CV boot")
            if (sexo == "Ambos") 
                pdf(file = paste(especie, "_cretal_plots.pdf", sep = ""), width = 9, 
                  height = 13)
            if (sexo == "Machos") 
                pdf(file = paste(especie, "_cretalmachos_plots.pdf", sep = ""), width = 9, 
                  height = 13)
            if (sexo == "Hembras") 
                pdf(file = paste(especie, "_cretalhembras_plots.pdf", sep = ""), width = 9, 
                  height = 13)
            layout(matrix(c(1, 1, 2, 3, 4, 0), 3, 2, byrow = T), 
                widths = c(3, 3), heights = c(1.5, 1.5))
            par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
                cex.lab = 1.2)
            plot(cre.tal.dat$tal ~ cre.tal.dat$eda, xlab = "Edad (años)", 
                ylab = paste("Talla (", unid, ")",sep=""), main = "Crecimiento en Talla", 
                pch = 1)
            lines(cre.tal.dat$eda, Li * (1 - exp(-cre.tal.fit[1] * 
                (cre.tal.dat$eda - cre.tal.fit[2]))), type = "l", 
                col = "red", lty = 2, lwd = 2)
            lines(cre.tal.dat$eda, Li * (1 - exp(-median(K.tal.b) * 
                (cre.tal.dat$eda - median(t0.b)))), type = "l", 
                col = "blue", lty = 3, lwd = 2)
            legend("bottomright", c("Predicción Original", 
                "Predicción Bootstrap"), lwd = 1, lty = 1, col = c("red", 
                "blue"),cex=0.8, bty="n")
            res.acf <- residuals(cre.tal.ori)
            x.acf <- sample(res.acf, length(res.acf), replace = F)
            acf(x.acf, main = "Correlograma")
            plot.ib(cre.tal.ori, which = 1, main = "Resíduos vs Predicciones")
            plot.ib(cre.tal.ori, which = 2, main = "Normal QQ Plot")
            mtext(outer = T, paste("Crecimiento en Talla: Gráficos Diagnóstico \nLinf (fija)=" , 
                Li, unid), side = 3, cex = 1.2)
            
            layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = F), widths = c(3, 
                3), heights = c(1.5, 1.5))
            par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
                cex.lab = 1.2)
            hist(K.tal.b, main = "BOOTSTRAP", xlab = "K (año-1)", 
                ylab = "Frecuencia")
            hist(t0.b, main = "BOOTSTRAP", xlab = "t0 (año)", 
                ylab = "Frecuencia")
            qqnorm(K.tal.b, xlab = "Cuantiles teóricos", ylab = "Cuantiles-K ", 
                main = "BOOTSTRAP")
            qqline(K.tal.b)
            qqnorm(t0.b, xlab = "Cuantiles teóricos", ylab = "Cuantiles-t0 ", 
                main = "BOOTSTRAP")
            qqline(t0.b)
            mtext(outer = T, paste("Crecimiento en Talla: Gráficos Diagnóstico\nLinf (fija)=" , 
                Li, unid), side = 3, cex = 1.2)
            dev.off()
        }
        if (sexo == "Ambos") 
            sink(file = paste(especie, "_cretal_resultados.txt", sep = ""))
        if (sexo == "Machos") 
            sink(file = paste(especie, "_cretalmachos_resultados.txt", sep = ""))
        if (sexo == "Hembras") 
            sink(file = paste(especie, "_cretalhembras_resultados.txt", sep = ""))
        cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", "Hora:", 
            format(Sys.time(), "%X"))
        cat("\n")
        cat("\n")
        cat("CRECIMIENTO EN TALLA: von Bertalanffy", "\n")
        cat("\n")
        cat("------------------------------------------------------", 
            "\n")
        cat("Especie:", especie, "\n")
        cat("Sexo:", sexo, "\n")
        cat("N:", ntal, "\n")
        if (Lfija == T) 
            cat(paste("Linf(fija) =", Li, unid), "\n")
        cat("\n")
        cat("Argumentos Rutina:", "\n")
        print(call)
        cat("------------------------------------------------------", 
            "\n")
        cat("\n")
        cat("RESULTADOS AJUSTE (datos originales)", "\n")
        cat("\n")
        print(cre.tal.ori,justify="right")
        cat("\n")
        cat("\n")
        cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP", "\n")
        cat("\n")
        print(par.tal,justify="right")
        sink()
    }
    cre.pes.f <- function(data, sexo) {
        cre.pes.dat <- data.frame(cbind(data$pes, data$eda))
        cre.pes.dat <- na.omit(cre.pes.dat)
        cre.pes.dat <- as.data.frame(cre.pes.dat)
        npes <- nrow(cre.pes.dat)
        dimnames(cre.pes.dat)[[2]] <- c("pes", "eda")
        if (Wfijo == F) {
          attach(cre.pes.dat)
          cre.pes.ori <- nls(pes ~ W * (1 - exp(-K * (eda - T0)))^b, start = list(W = Wi, K = Kwi, T0 = T0wi), 
                               control = nls.control(maxiter = 100, tol = 1e-04))
            detach(cre.pes.dat)
            cre.pes.fit <- coef(cre.pes.ori)
            cre.pes.fit <- as.vector(cre.pes.fit)
            names(cre.pes.fit) <- c("Winf (g)", "k (año-1)", 
                "t0 (año)")
            Wi.ori <- cre.pes.fit[1]
            K.pes.ori <- cre.pes.fit[2]
            t0.pes.ori <- cre.pes.fit[3]
            cre.pes.fun <- function(data) coef(nls(pes ~ 
                W * (1 - exp(-K * (eda - T0)))^b, start = coef(cre.pes.ori), 
                control = nls.control(maxiter = 100, tol = 1e-04), data=data))
            cre.pes.case <- function(data, i) cre.pes.fun(data[i,])
            cre.pes.boot <- boot(cre.pes.dat, cre.pes.case, R = n)
            Wi.b <- cre.pes.boot$t[, 1]
            K.pes.b <- cre.pes.boot$t[, 2]
            t0.pes.b <- cre.pes.boot$t[, 3]
            Wi.boot <- as.vector(cv(Wi.b))
            K.pes.boot <- as.vector(cv(K.pes.b))
            t0.pes.boot <- as.vector(cv(t0.pes.b))
            par.pes.boot <- rbind(Wi.boot, K.pes.boot, t0.pes.boot)
            par.pes.ori <- round(coef(cre.pes.ori), 6)
            par.pes <- cbind(par.pes.ori, par.pes.boot)
            dimnames(par.pes)[[1]] <- c("Winf(g)", "k(año-1)", "t0(año)")
            dimnames(par.pes)[[2]] <- c("Estima original", "Estima boot", 
                "CV boot")
            if (sexo == "Ambos") 
                pdf(file = paste(especie, "_crepeso_plots.pdf", sep = ""), width = 9, 
                  height = 13)
            if (sexo == "Machos") 
                pdf(file = paste(especie, "_crepesomachos_plots.pdf", sep = ""), width = 9, 
                  height = 13)
            if (sexo == "Hembras") 
                pdf(file = paste(especie, "_crepesohembras_plots.pdf", sep = ""), 
                  width = 9, height = 13)
            layout(matrix(c(1, 1, 2, 3, 4, 0), 3, 2, byrow = T), 
                widths = c(1, 1), heights = c(1, 1))
            par(oma = c(1, 1, 2, 1), font = 3, cex = 1, cex.main = 1.2, 
                cex.lab = 1.2)
            plot(cre.pes.dat$pes ~ cre.pes.dat$eda, xlab = "Edad (años)", 
                ylab = "Peso (g)", main = "Crecimiento en Peso", 
                pch = 1)
            lines(cre.pes.dat$eda, cre.pes.fit[1] * (1 - exp(-cre.pes.fit[2] * 
                (cre.pes.dat$eda - cre.pes.fit[3])))^b, type = "l", 
                col = "red", lty = 2, lwd = 2)
            lines(cre.pes.dat$eda, median(Wi.b) * (1 - exp(-median(K.pes.b) * 
                (cre.pes.dat$eda - median(t0.pes.b))))^b, type = "l", 
                col = "blue", lty = 3, lwd = 2)
            legend("bottomright", c("Predicción Original", 
                "Predicción Bootstrap"), lwd = 1, lty = 1, col = c("red", 
                "blue"),cex=0.8,bty="n")
            res.acf <- residuals(cre.pes.ori)
            x.acf <- sample(res.acf, length(res.acf), replace = F)
            acf(x.acf, main = "Correlograma")
            plot.ib(cre.pes.ori, which = 1, main = "Resíduos vs Predicciones")
            plot.ib(cre.pes.ori, which = 2, main = "Normal QQ Plot")
            mtext(outer = T, " Crecimiento en Peso: Gráficos Diagnóstico" , 
                side = 3, cex = 1.2)
            layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = F), 
                widths = c(3, 3), heights = c(1.5, 1.5))
            par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
                cex.lab = 1.2)
            hist(Wi.b, main = "BOOTSTRAP", xlab = "Winf (g)", 
                ylab = "Frecuencia")
            hist(K.pes.b, main = "BOOTSTRAP", xlab = "K (año-1)", 
                ylab = "Frecuencia")
            hist(t0.pes.b, main = "BOOTSTRAP", xlab = "t0 (año)", 
                ylab = "Frecuencia")
            qqnorm(Wi.b, xlab = "Cuantiles teóricos", ylab = "Cuantiles-Winf ", 
                main = "BOOTSTRAP")
            qqline(Wi.b)
            qqnorm(K.pes.b, xlab = "Cuantiles teóricos", ylab = "Cuantiles-K ", 
                main = "BOOTSTRAP")
            qqline(K.pes.b)
            qqnorm(t0.pes.b, xlab = "Cuantiles teóricos", ylab = "Cuantiles-t0 ", 
                main = "BOOTSTRAP")
            qqline(t0.pes.b)
            mtext(outer = T, " Crecimiento en Peso: Gráficos Diagnóstico" , 
                side = 3, cex = 1.2)
            dev.off()
        }
        if (Wfijo == T) {
            attach(cre.pes.dat)
            cre.pes.ori <- nls(pes ~ Wi * (1 - exp(-K * (eda - 
                T0)))^b, start = list(K = Kwi, T0 = T0wi), control = nls.control(maxiter = 100, 
                tol = 1e-04))
            detach(cre.pes.dat)
            cre.pes.fit <- coef(cre.pes.ori)
            cre.pes.fit <- as.vector(cre.pes.fit)
            names(cre.pes.fit) <- c("K (año-1)", "t0 (año)")
            K.pes.ori <- cre.pes.fit[1]
            t0.pes.ori <- cre.pes.fit[2]
            cre.pes.fun <- function(data) coef(nls(pes ~ 
                Wi * (1 - exp(-K * (eda - T0)))^b, start = coef(cre.pes.ori), 
                control = nls.control(maxiter = 100, tol = 1e-04), data=data))
            cre.pes.case <- function(data, i) cre.pes.fun(data[i, ])
            cre.pes.boot <- boot(cre.pes.dat, cre.pes.case, R = n)
            K.pes.b <- cre.pes.boot$t[, 1]
            t0.pes.b <- cre.pes.boot$t[, 2]
            K.pes.boot <- as.vector(cv(K.pes.b))
            t0.pes.boot <- as.vector(cv(t0.pes.b))
            par.pes.boot <- rbind(K.pes.boot, t0.pes.boot)
            par.pes.ori <- round(coef(cre.pes.ori), 6)
            par.pes <- cbind(par.pes.ori, par.pes.boot)
            dimnames(par.pes)[[1]] <- c("k(año-1)", "t0(año)")
            dimnames(par.pes)[[2]] <- c("Estima original", "Estima boot", 
                "CV boot")
            if (sexo == "Ambos") 
                pdf(file = paste(especie, "_crepeso_plots.pdf", sep = ""), width = 9, 
                  height = 13)
            if (sexo == "Machos") 
                pdf(file = paste(especie, "_crepesomachos_plots.pdf", sep = ""), width = 9, 
                  height = 13)
            if (sexo == "Hembras") 
                pdf(file = paste(especie, "_crepesohembras_plots.pdf", sep = ""), 
                  width = 9, height = 13)
            layout(matrix(c(1, 1, 2, 3, 4, 0), 3, 2, byrow = T), 
                widths = c(1, 1), heights = c(1, 1))
            par(oma = c(1, 1, 2, 1), font = 3, cex = 1, cex.main = 1.2, 
                cex.lab = 1.2)
            plot(cre.pes.dat$pes ~ cre.pes.dat$eda, xlab = "Edad (años)", 
                ylab = "Peso (g)", main = "Crecimiento en Peso", 
                pch = 1)
            lines(cre.pes.dat$eda, Wi * (1 - exp(-cre.pes.fit[1] * 
                (cre.pes.dat$eda - cre.pes.fit[2])))^b, type = "l", 
                col = "red", lty = 2, lwd = 2)
            lines(cre.pes.dat$eda, Wi * (1 - exp(-median(K.pes.b) * 
                (cre.pes.dat$eda - median(t0.pes.b))))^b, type = "l", 
                col = "blue", lty = 3, lwd = 2)
            legend("bottomright",c("Predicción Original", 
                "Predicción Bootstrap"), lwd = 1, lty = 1, col = c("red", 
                "blue"), cex=0.8, bty="n")
            res.acf <- residuals(cre.pes.ori)
            x.acf <- sample(res.acf, length(res.acf), replace = F)
            acf(x.acf, main = "Correlograma")
            plot.ib(cre.pes.ori, which = 1, main = "Resíduos vs Predicciones")
            plot.ib(cre.pes.ori, which = 2, main = "Normal QQ Plot")
            mtext(outer = T, paste("Crecimiento en Peso: Gráficos Diagnóstico\nWinf (fijo)=" , 
                Wi, "g"), side = 3, cex = 1.2)
            layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = F), widths = c(3, 
                3), heights = c(1.5, 1.5))
            par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
                cex.lab = 1.2)
            hist(K.pes.b, main = "BOOTSTRAP", xlab = "K (año-1)", 
                ylab = "Frecuencia")
            hist(t0.pes.b, main = "BOOTSTRAP", xlab = "t0 (año)", 
                ylab = "Frecuencia")
            qqnorm(K.pes.b, xlab = "Cuantiles teóricos", ylab = "Cuantiles-K ", 
                main = "BOOTSTRAP")
            qqline(K.pes.b)
            qqnorm(t0.pes.b, xlab = "Cuantiles teóricos", ylab = "Cuantiles-t0 ", 
                main = "BOOTSTRAP")
            qqline(t0.pes.b)
            mtext(outer = T, paste("Crecimiento en Peso: Gráficos Diagnóstico\nWinf (fijo)=" , 
                Wi, "g"), side = 3, cex = 1.2)
            dev.off()
        }
        if (sexo == "Ambos") 
            sink(file = paste(especie, "_crepeso_resultados.txt", sep = ""))
        if (sexo == "Machos") 
            sink(file = paste(especie, "_crepesomachos_resultados.txt", sep = ""))
        if (sexo == "Hembras") 
            sink(file = paste(especie, "_crepesohembras_resultados.txt", sep = ""))
        cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", "Hora:", 
            format(Sys.time(), "%X"))
        cat("\n")
        cat("\n")
        cat("CRECIMIENTO EN PESO: von Bertalanffy", "\n")
        cat("\n")
        cat("------------------------------------------------------------", 
            "\n")
        cat("Especie:", especie, "\n")
        cat("Sexo:", sexo, "\n")
        cat("N:", npes, "\n")
        if (Wfijo == T) 
            cat(paste("Winf(fijo) =", Wi, "g"), "\n")
        cat("\n")
        cat("Argumentos Rutina:", "\n")
        print(call)
        cat("------------------------------------------------------------", 
            "\n")
        cat("\n")
        cat("RESULTADOS AJUSTE (datos originales)", "\n")
        cat("\n")
        print(cre.pes.ori,justify="right")
        cat("\n")
        cat("\n")
        cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP", "\n")
        cat("\n")
        print(par.pes, justify="right")
        sink()
    }
    dat <- read.csv(file = paste(especie,".csv", sep = ""), header = T)
    dat <- data.frame(dat)
    if (cl == 0.5) 
        dat$tal <- mcinf(dat$tal)
    if (cl == 1) 
        dat$tal <- floor(dat$tal)
    dat$tal <- dat$tal + (cl/2)
    dat$eda <- dat$eda + 0.5
    orde <- sort.list(dat[, 3])
    dat <- dat[orde, ]
    if (sex == F) {
        cre.tal.f(dat, "Ambos")
        cre.pes.f(dat, "Ambos")
    }
    if (sex == T) {
        datm <- dat[dat[, 5] == 1, ]
        datf <- dat[dat[, 5] == 2, ]
        rm(dat)
        cre.tal.f(datm, "Machos")
        cre.tal.f(datf, "Hembras")
        cre.pes.f(datm, "Machos")
        cre.pes.f(datf, "Hembras")
    }
    path <- getwd()
    cat("\n")
    cat(paste("EL PROGRAMA HA TERMINADO.\n","PUEDES VER LOS RESULTADOS EN EL DIRECTORIO DE TRABAJO:", path))
    cat("\n")
}

