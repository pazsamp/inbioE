"madurez" <-
function (especie = "nombre especie", cl = 1, unid = "cm", edad = T, 
    sex = F, n = 1000) 
{
    call <- match.call()
    cv <- function(data) {
        data <- as.numeric(data)
        med <- round(median(data), 4)
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
    ilogit <- function(l) {
        exp(l)/(1 + exp(l))
    }
    mad.tal <- function(data, sexo) {
        mad.tal.dat <- data.frame(cbind(data$tal, data$mad))
        mad.tal.dat <- na.omit(mad.tal.dat)
        mad.tal.dat <- as.data.frame(mad.tal.dat)
        ntal <- nrow(mad.tal.dat)
        dimnames(mad.tal.dat)[[2]] <- c("tal", "mad")
        attach(mad.tal.dat)
        p.t <- table(as.factor(tal), mad)
        p2.t <- p.t[, 2]/(p.t[, 1] + p.t[, 2])
        x.p2.t <- as.numeric(levels(as.factor(tal)))
        mad.tal.ori <- glm(mad ~ tal, family = binomial(link = logit), 
            data = mad.tal.dat)
        mad.tal.fun <- function(data) coef(glm(mad ~ tal, family = binomial(link = logit), 
            data = data))
        mad.tal.case <- function(data, i) mad.tal.fun(data[i, 
            ])
        mad.tal.boot <- boot(mad.tal.dat, mad.tal.case, R = n)
        b0 <- mad.tal.boot$t[, 1]
        b1 <- mad.tal.boot$t[, 2]
        b0.boot <- cv(b0)
        b1.boot <- cv(b1)
        L50 <- (-b0/b1)
        L50.boot <- cv(L50)
        par.boot <- rbind(b0.boot, b1.boot, L50.boot)
        par.ori <- round(c(coef(mad.tal.ori), -coef(mad.tal.ori)[1]/coef(mad.tal.ori)[2]), 
            4)
        par <- cbind(par.ori, par.boot)
        dimnames(par)[[1]] <- c("B0", "B1", paste("L50(", unid, ")",sep=""))
        dimnames(par)[[2]] <- c("Estima datos", "Estima boot", 
            "CV boot")
        detach(mad.tal.dat)
        if (sexo == "A") 
            pdf(file = paste(especie, "_madtal_plots.pdf", sep = ""), width = 9, 
                height = 13)
        if (sexo == "M") 
            pdf(file = paste(especie, "_madtalmachos_plots.pdf", sep = ""), width = 9, 
                height = 13)
        if (sexo == "F") 
        pdf(file = paste(especie, "_madtalhembras_plots.pdf", sep = ""), width = 9, 
                height = 13)
        layout(matrix(c(1, 1, 2, 3, 4, 5), 3, 2, byrow = T), 
            widths = c(3, 3), heights = c(1.5, 1.5))
        par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
            cex.lab = 1.2)
        plot(p2.t ~ x.p2.t, xlab = paste("Talla (", unid, ")", sep=""), 
            ylab = " Proporción     Madurez", main = "Curva de Madurez")
        rangt <- seq(min(mad.tal.dat$tal), max(mad.tal.dat$tal), 
            by = 0.2)
        lines(rangt, ilogit(mad.tal.ori$coef[1] + rangt * mad.tal.ori$coef[2]), 
            col = "red", lty = 2, lwd = 2)
        lines(rangt, ilogit(median(b0) + rangt * median(b1)), 
            col = "blue", lty = 3, lwd = 2)
        x0.t <- (-mad.tal.ori$coef[1])/mad.tal.ori$coef[2]
        segments(-0.5, 0.5, x0.t, 0.5, lty = 3)
        segments(x0.t, 0.5, x0.t, -1, lty = 3)
        legend("bottomright", c("Predicción Original", "Predicción Bootstrap"), 
        col = c("red", "blue"), lty = c(2, 3), lwd = c(2, 
        2), bty="n")
        
        res.acf <- residuals(mad.tal.ori)
        x.acf <- sample(res.acf, length(res.acf), replace = F)
        acf(x.acf, main = "Correlograma")
        r.t <- ilogit(mad.tal.ori$coef[1] + x.p2.t * mad.tal.ori$coef[2])
        res.t <- (r.t - p2.t)
        res.t <- (res.t - mean(res.t))/sd(res.t)
        plot(res.t ~ r.t, ylab = "Resíduos estandarizados", xlab = "Predicciones", 
            main = "Original: Análisis de Resíduos" )
        abline(h = 0, lty = 3)
        hist(L50, main = "BOOTSTRAP", xlab = paste("L50 (", unid, 
            ")", sep=""), ylab = "Frecuencia")
        qqnorm(L50, xlab = "Cuantiles teóricos", ylab = "Cuantiles-L50 ", 
            main = "BOOTSTRAP")
        qqline(L50)
        mtext(outer = T, "Talla de Madurez: Gráficos Diagnóstico" , 
            side = 3, cex = 1.2)
        graphics.off()
        res <- list(ntal = ntal, glm.ori = mad.tal.ori, b0.boot = b0.boot, 
            b1.boot = b1.boot, L50.boot = L50.boot, parametro = par)
        rm(mad.tal.dat, mad.tal.ori, mad.tal.boot, mad.tal.case, 
            mad.tal.fun, b0, b1, L50)
        res
    }
    mad.eda <- function(data, sexo) {
        mad.eda.dat <- data.frame(cbind(data$eda, data$mad))
        mad.eda.dat <- na.omit(mad.eda.dat)
        mad.eda.dat <- as.data.frame(mad.eda.dat)
        neda <- nrow(mad.eda.dat)
        dimnames(mad.eda.dat)[[2]] <- c("eda", "mad")
        attach(mad.eda.dat)
        p.e <- table(as.factor(eda), mad)
        p2.e <- p.e[, 2]/(p.e[, 1] + p.e[, 2])
        x.p2.e <- as.numeric(levels(as.factor(eda)))
        mad.eda.ori <- glm(mad ~ eda, family = binomial(link = logit), 
            data = mad.eda.dat)
        mad.eda.fun <- function(data) coef(glm(mad ~ eda, family = binomial(link = logit), 
            data = data))
        mad.eda.case <- function(data, i) mad.eda.fun(data[i, 
            ])
        mad.eda.boot <- boot(mad.eda.dat, mad.eda.case, R = n)
        b0 <- mad.eda.boot$t[, 1]
        b1 <- mad.eda.boot$t[, 2]
        b0.boot <- cv(b0)
        b1.boot <- cv(b1)
        E50 <- (-b0/b1)
        E50.boot <- cv(E50)
        detach(mad.eda.dat)
        par.boot <- rbind(b0.boot, b1.boot, E50.boot)
        par.ori <- round(c(coef(mad.eda.ori), -coef(mad.eda.ori)[1]/coef(mad.eda.ori)[2]), 
            4)
        par <- cbind(par.ori, par.boot)
        dimnames(par)[[1]] <- c("B0", "B1", "E50(año)")
        dimnames(par)[[2]] <- c("Estima datos", "Estima boot", "CV boot")
        if (sexo == "A") 
            pdf(file = paste(especie, "_madedad_plots.pdf", sep = ""), width = 9, 
                height = 13)
        if (sexo == "M") 
            pdf(file = paste(especie, "-madedadmachos_plots.pdf", sep = ""), width = 9, 
                height = 13)
        if (sexo == "F") 
            pdf(file = paste(especie, "_madedadhembras_plots.pdf", sep = ""), width = 9, 
                height = 13)
        layout(matrix(c(1, 1, 2, 3, 4, 5), 3, 2, byrow = T), 
            widths = c(3, 3), heights = c(1.5, 1.5))
        par(oma = c(2, 2, 3, 2), font = 3, cex = 1, cex.main = 1.2, 
            cex.lab = 1.2)
        plot(p2.e ~ x.p2.e, xlab = "Edad (año)", ylab = " Proporción  Madurez", 
            main = "Curva de  Madurez")
        rang <- seq(min(mad.eda.dat$eda), max(mad.eda.dat$eda), 
            by = 0.2)
        lines(rang, ilogit(mad.eda.ori$coef[1] + rang * mad.eda.ori$coef[2]), 
            col = "red", lty = 2, lwd = 2)
        lines(rang, ilogit(median(b0) + rang * median(b1)), col = "blue", 
            lty = 3, lwd = 2)
        x0.e <- (-mad.eda.ori$coef[1])/mad.eda.ori$coef[2]
        segments(-0.5, 0.5, x0.e, 0.5, lty = 3)
        segments(x0.e, 0.5, x0.e, -1, lty = 3)
        legend("bottomright", c("Predicción Original", "Predicción Bootstrap"), 
            col = c("red", "blue"), lty = c(2, 3), lwd = c(2, 
                2), bty="n")
        res.eda <- residuals(mad.eda.ori)
        x.acf.eda <- sample(res.eda, length(res.eda), replace = F)
        acf(x.acf.eda, main = "Correlograma")
        r.e <- ilogit(mad.eda.ori$coef[1] + x.p2.e * mad.eda.ori$coef[2])
        res.e <- (r.e - p2.e)
        res.e <- (res.e - mean(res.e))/sd(res.e)
        plot(res.e ~ r.e, ylab = "Resíduos estandarizados", xlab = " Predicciones ", 
            main = "Original: Análisis de Resíduos" )
        abline(h = 0, lty = 3)
        hist(E50, main = "BOOTSTRAP", xlab = "E50(año)", ylab = "Frecuencia")
        qqnorm(E50, xlab = "Cuantiles teóricos", ylab = "Cuantiles-E50 ", 
            main = "BOOTSTRAP")
        qqline(E50)
        mtext(outer = T, "Edad de Madurez: Gráficos Diagnóstico" , 
            side = 3, cex = 1.2)
        graphics.off()
        res <- list(neda = neda, glm.ori = mad.eda.ori, b0.boot = b0.boot, 
            b1.boot = b1.boot, E50.boot = E50.boot, parametro = par)
        rm(mad.eda.dat, mad.eda.ori, mad.eda.boot, mad.eda.case, 
            mad.eda.fun, b0, b1, E50)
        res
    }
    dat <- read.csv(file = paste(especie, ".csv", sep = ""), sep = ",", header = T)
    dat <- data.frame(dat)
    if (cl == 0.5) 
        dat$tal <- mcinf(dat$tal)
    if (cl == 1) 
        dat$tal <- floor(dat$tal)
    dat$tal <- dat$tal + (cl/2)
    if (edad == T) 
        dat$eda <- dat$eda + 0.5
    if (sex == F) {
        mad.tal.res <- mad.tal(dat, "A")
        sink(file = paste(especie,"_madtal_resultados.txt", sep = ""))
        cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", "Hora:", 
            format(Sys.time(), "%X"))
        cat("\n")
        cat("\n")
        cat("\n")
        cat("TALLA DE MADUREZ", "\n")
        cat("\n")
        cat("-------------------------------------------------------", 
            "\n")
        cat("Especie:", especie, "\n")
        cat("Sexo: Ambos", "\n")
        cat("N:", mad.tal.res$ntal, "\n")
        cat("\n")
        cat("Argumentos Rutina:", "\n")
        print(call)
        cat("-------------------------------------------------------", 
            "\n")
        cat("\n")
        cat("RESULTADOS AJUSTE(datos originales)", "\n")
        cat("\n")
        print(summary(mad.tal.res$glm.ori))
        cat("\n")
        cat("\n")
        cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP", "\n")
        cat("\n")
        print(mad.tal.res$parametro)
        sink()
        rm(mad.tal.res)
        if (edad == T) {
            mad.eda.res <- mad.eda(dat, "A")
            sink(file = paste(especie, "_madedad_resultados.txt", sep = ""))
            cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", 
                "Hora:", format(Sys.time(), "%X"))
            cat("\n")
            cat("\n")
            cat("\n")
            cat("EDAD DE MADUREZ", "\n")
            cat("\n")
            cat("-------------------------------------------------------", 
                "\n")
            cat("Especie:", especie, "\n")
            cat("Sexo: Ambos", "\n")
            cat("N:", mad.eda.res$neda, "\n")
            cat("\n")
            cat("Argumentos Rutina:", "\n")
            print(call)
            cat("-------------------------------------------------------", 
                "\n")
            cat("\n")
            cat("RESULTADOS AJUSTE(datos originales)", "\n")
            cat("\n")
            print(summary(mad.eda.res$glm.ori))
            cat("\n")
            cat("\n")
            cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP", 
                "\n")
            cat("\n")
            print(mad.eda.res$parametro)
            sink()
            rm(mad.eda.res)
        }
    }
    if (sex == T) {
        datm <- dat[dat$sex == 1, ]
        datf <- dat[dat$sex == 2, ]
        rm(dat)
        mad.tal.resf <- mad.tal(datf, "F")
        mad.tal.resm <- mad.tal(datm, "M")
        sink(file = paste(especie, "_madtalhembras_resultados.txt", sep = ""))
        cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", "Hora:", 
            format(Sys.time(), "%X"))
        cat("\n")
        cat("\n")
        cat("\n")
        cat("TALLA DE MADUREZ", "\n")
        cat("\n")
        cat("-------------------------------------------------------", 
            "\n")
        cat("Especie:", especie, "\n")
        cat("Sexo: Hembras", "\n")
        cat("N:", mad.tal.resf$ntal, "\n")
        cat("\n")
        cat("Argumentos Rutina:", "\n")
        print(call)
        cat("-------------------------------------------------------", 
            "\n")
        cat("\n")
        cat("RESULTADOS AJUSTE (datos originales)", "\n")
        print(summary(mad.tal.resf$glm.ori))
        cat("\n")
        cat("\n")
        cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP", "\n")
        cat("\n")
        print(mad.tal.resf$parametro)
        sink()
        sink(file = paste(especie,"_madtalmachos_resultados.txt", sep = ""))
        cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", "Hora:", 
            format(Sys.time(), "%X"))
        cat("\n")
        cat("\n")
        cat("\n")
        cat("TALLA DE MADUREZ", "\n")
        cat("\n")
        cat("-------------------------------------------------------", 
            "\n")
        cat("Especie:", especie, "\n")
        cat("Sexo: Machos", "\n")
        cat("N:", mad.tal.resm$ntal, "\n")
        cat("\n")
        cat("Argumentos Rutina:", "\n")
        print(call)
        cat("-------------------------------------------------------", 
            "\n")
        cat("\n")
        cat("RESULTADOS AJUSTE(datos originales)", "\n")
        print(summary(mad.tal.resm$glm.ori))
        cat("\n")
        cat("\n")
        cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP", "\n")
        cat("\n")
        print(mad.tal.resm$parametro)
        sink()
        rm(mad.tal.resf, mad.tal.resm)
        if (edad == T) {
            mad.eda.resf <- mad.eda(datf, "F")
            mad.eda.resm <- mad.eda(datm, "M")
            sink(file = paste(especie, "_madedadhembras_resultados.txt", sep = ""))
            cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", 
                "Hora:", format(Sys.time(), "%X"))
            cat("\n")
            cat("\n")
            cat("\n")
            cat("EDAD DE MADUREZ", "\n")
            cat("\n")
            cat("-------------------------------------------------------", 
                "\n")
            cat("Especie:", especie, "\n")
            cat("Sexo: Hembras", "\n")
            cat("N:", mad.eda.resf$neda, "\n")
            cat("\n")
            cat("Argumentos Rutina:", "\n")
            print(call)
            cat("-------------------------------------------------------", 
                "\n")
            cat("\n")
            cat("RESULTADOS AJUSTE(datos originales)", "\n")
            cat("\n")
            print(summary(mad.eda.resf$glm.ori))
            cat("\n")
            cat("\n")
            cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP", 
                "\n")
            cat("\n")
            print(mad.eda.resf$parametro)
            sink()
            sink(file = paste(especie, "_madedadmachos_resultados.txt", sep = ""))
            cat("Día:", format(Sys.time(), "%d-%b-%Y"), "   ", 
                "Hora:", format(Sys.time(), "%X"))
            cat("\n")
            cat("\n")
            cat("\n")
            cat("EDAD DE MADUREZ", "\n")
            cat("\n")
            cat("-------------------------------------------------------", 
                "\n")
            cat("Especie:", especie, "\n")
            cat("Sexo: Machos", "\n")
            cat("N:", mad.eda.resm$neda, "\n")
            cat("\n")
            cat("Argumentos Rutina:", "\n")
            print(call)
            cat("-------------------------------------------------------", 
                "\n")
            cat("\n")
            cat("RESULTADOS AJUSTE (datos originales)", "\n")
            cat("\n")
            print(summary(mad.eda.resm$glm.ori))
            cat("\n")
            cat("\n")
            cat("PARÁMETROS ESTIMADOS: ORIGINALES-BOOTSTRAP", 
                "\n")
            cat("\n")
            print(mad.eda.resm$parametro)
            sink()
            rm(mad.eda.resf, mad.eda.resm)
        }
    }
    
    path <- getwd()
    cat("\n")
    cat(paste("EL PROGRAMA HA TERMINADO.\n","PUEDES VER LOS RESULTADOS EN EL DIRECTORIO DE TRABAJO:", path))
    cat("\n")
}

