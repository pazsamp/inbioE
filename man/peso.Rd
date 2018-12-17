\name{peso}
\alias{peso}
\title{Relacion Talla-Peso}
\description{
Calcula la relación talla-peso por estimación no linear y ajuste de mínimos cuadrados. 
Mediante bootstrap no paramétrico estima el coeficiente de variación de los parámetros de la relación talla-peso.
}
\usage{
peso(especie = "nombre especie", cl = 1, unid = "cm", sex = F, a = 5e-04, b = 3, n = 1000)
}
\arguments{
  \item{especie}{Nombre del archivo donde están los datos de entrada, sin la extensión y sin el directorio.}
  \item{cl}{Intervalo de talla. El valor por defecto es 1.}
  \item{unid}{Unidad de los datos de talla: cm o mm. Por defecto unid="cm".}
  \item{sex}{Si sex=F el proceso analiza los datos de machos y hembras conjuntamente. Si sex=T el proceso analiza separadamente machos y hembras. Por defecto sex=F.}
  \item{a}{Valor de partida ("start value" o "semilla") de la ordenada.}
  \item{b}{Valor de partida ("start value" o "semilla") de la pendiente.}
  \item{n}{Número de réplicas bootstrap. El número recomendado es 1000.}
}

\details{
   Es importante que los valores de partida ("start value" o "semilla") no sean muy diferentes de los que se quiere calcular (se pueden consultar valores de otros años o zonas), si no, el proceso
   iterativo por el que se realiza el ajuste podría fallar.  
   }

\value{
Archivo "txt" con los resultados de los datos originales y de las réplicas bootstrap, así como su coeficiente de variación.

Archivo "pdf" con las gráficas del ajuste para los datos originales y para las réplicas bootstrap, gráficas del análisis de resíduos del ajuste de los datos originales y gráficas de análisis de las réplicas bootstrap.

Si el proceso se realiza separadamente para machos y hembras, los archivos de salida también son independientes.

}
\author{Paz Sampedro}

\references{Efron, B. 1979. Bootstrap methods: another look at the jackknife. Annals of Statistics, 7. 1-26 p.}

\examples{
## Requiere que exista el archivo "merluza.csv" en el directorio de trabajo
## Calcula la relación talla-peso separadamente para machos y hembras
   peso (especie="merluza", sex=T, a= 0.005, b=3)

## Requiere que exista el archivo "choco.csv" en el directorio de trabajo
## Calcula la relación talla-peso para una especie con intervalo de talla de 0.5 cm, separadamente por sexos
   peso(especie="choco", cl=0.5, sex=T, a=0.05, b=2.4)
 } 
