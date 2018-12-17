\name{sexratio}
\alias{sexratio}
\title{Sexratio por Talla y Edad}
\description{
Calcula el porcentaje de hembras por talla y por edad, y su coeficiente de variación mediante bootstrap no paramétrico.
** Sexratio= (hembras/(hembras+machos)
}
\usage{
sexratio(especie = "nombre especie", cl = 1, unid = "cm", edad = T, n = 1000)
}

\arguments{
  \item{especie}{Nombre del archivo donde están los datos de entrada, sin la extensión y sin el directorio.}
  \item{cl}{Intervalo de talla usado en la sexratio por talla. El valor por defecto es 1.}
  \item{unid}{Unidad de los datos de talla: cm o mm. Por defecto unid="cm".}
  \item{edad}{Si edad=T la rutina calcula el sexratio por edad; Si edad=F la rutina no calcula sexratio por edad. Por defecto edad = T.}
  \item{n}{Número de réplicas bootstrap. El número recomendado es 1000.}
}

\value{
Archivo "txt" con los resultados de sexratio por talla y su coeficiente de variación global y los resultados para sexratio por edad y coeficiente de variación global.

Archivo "pdf" con dos gráficas, una gráfico del sexratio y coeficiente de variación por talla y otra del sexratio y coeficiente de variación por edad.
}
 
\author{Paz Sampedro}

\references{Efron, B. 1979. Bootstrap methods: another look at the jackknife. Annals of Statistics, 7. 1-26 p.}

\examples{
# Requiere que exista el archivo "merluza.csv" en el directorio de trabajo
# Calcula el sexratio por talla y el sexratio por edad
  sexratio(especie="merluza", edad=T)

# Requiere que exista el archivo "cigala.csv" en el directorio de trabajo
# Calcula el sexratio por talla 
  sexratio(especie="cigala", edad=F)
}
