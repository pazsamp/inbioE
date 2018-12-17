\name{madurez}
\alias{madurez}
\title{Madurez Sexual en Talla y Edad}
\description{
Calcula la talla y la edad de madurez sexual aplicando un glm con errores binomiales. Mediante bootstrap no paramétrico estima los coeficientes de variación de: los parámetros de la función, la L50 y la A50.
}
\usage{
madurez(especie = "nombre especie", cl = 1, unid = "cm", edad = T, sex = F, n = 1000)
}
\arguments{
  \item{especie}{Nombre del archivo donde están los datos de entrada, sin la extensión y sin el directorio.}
  \item{cl}{Intervalo de talla de la especie. El valor por defecto es 1.}
  \item{unid}{Unidad de los datos de talla: cm o mm. Por defecto unid="cm"}
  \item{edad}{Si edad=T, se calcula la edad de madurez. Si edad=F, el proceso no calcula la edad de madurez. Por defecto edad=T.}
  \item{sex}{Si sex=F, el proceso analiza los datos de machos y hembras conjuntamente. Si sex=T, se analizan separadamente machos y hembras. Por defecto sex=F.}
  \item{n}{Número de réplicas bootstrap. El número recomendado es 1000.}
}

\value{
  Dos archivo "txt" con los resultados de los datos originales y de las réplicas bootstrap y su coeficiente de variación,uno para madurez en talla y otro para madurez en edad.

  Dos archivos "pdf" con las gráficas del ajuste con datos originales y bootstrap, las gráficas del análisis de resíduos del ajuste de los datos originales y análisis de las réplicas bootstrap, un archivo para madurez en talla y otro para madurez 
  en edad.           

  Si el proceso se realiza separadamente para machos y hembras, los archivos de salida también son independientes.

}
\references{
  Davison, A.C. & D.V. Hinkley. 1999. Bootstrap methods and their application. Cambridge University Press, 575 pp.

  Efron, B. 1979. Bootstrap methods: another look at the jackknife. Annals of Statistics, 7. 1-26 p.

  Roa, R., Ernst, B & F. Tapia. 1999. Estimation of size at sexual maturity: an evaluation of analytical and resampling procedures. Fishery Bulletin 97: 570-580.

}

\author{Paz Sampedro}

\examples{
## Requiere que exista el archivo "lepwhi2004.csv" en el directorio de trabajo
## Calcula la talla y la edad de madurez sexual y su coeficiente de variación conjuntamente para ambos sexos
   madurez (especie="lepwhi2004", sex=F)

## Requiere que exista el archivo "cangrejoreal.csv" en el directorio de trabajo
## Calcula sólo la talla de madurez sexual (edad=F) y su coeficiente de variación separadamente para machos y 
## hembras. La talla es en milímetros.
   madurez(especie="cangrejoreal", unid="mm", edad=F, sex=T)
}

