\name{crecimiento}
\alias{crecimiento}
\title{Crecimiento en Talla y en Peso}
\description{
Ajusta la ecuación de crecimiento en talla y en peso de von Bertalanffy. Mediante bootstrap no paramétrico estima el coeficiente de variación de los parámetros de la ecuación.
}
\usage{
crecimiento(especie = "nombre especie", cl = 1, unid = "cm", sex = F, b = 3, Li = 80, Lfija = F, Ki = 0.7, T0i = -1, Wi = 2500, Wfijo = F, Kwi = 0.4, T0wi = -1, n = 1000)
}
\arguments{
  \item{especie}{Archivo donde están los datos de entrada, sin la extensión y sin el directorio.}
  \item{cl}{Intervalo de talla de la especie.El valor por defecto es 1.}
  \item{unid}{Unidad de los datos de talla: cm o mm. Por defecto unid="cm".}
  \item{sex}{Si sex=F, el proceso analiza los datos de machos y hembras conjuntamente. Si sex=T, analiza separadamente machos y hembras. Por defecto sex=F.}
  \item{b}{Parámetro de la relación talla-peso empleado en el crecimiento en peso. Por defecto b=3.}
  \item{Li}{Valor de partida ("start value" o "semilla") de la Linf (longitud asintótica) que usa el proceso para iniciar el ajuste del crecimiento en talla.}
  \item{Lfija}{Si Lfija=T, fija el valor dado en "Li" como Linf en la ecuación de von Bertalanffy. Si Lfija=F, calcula el parámetro Linf y su coeficiente de variación. Por defecto Lfija=F.}
  \item{Ki}{Valor de partida ("start value" o "semilla") de la Ki (parámetro de curvatura) que usa el proceso para iniciar el ajuste del crecimiento en talla.}
  \item{T0i}{Valor de partida ("start value" o "semilla") de la T0i (punto en el tiempo, en el que el pez tiene una talla cero) que usa el proceso para iniciar el ajuste del crecimiento en talla.}
  \item{Wi}{Valor de partida de ("start value" o "semilla") del Winf (peso asintótico) que usa el proceso para iniciar el ajuste del crecimiento en peso.}
  \item{Wfijo}{Si Wfijo=T, fija el valor dado en "Wi" como Winf en la ecuación de von Bertalanffy. Si Wfijo=F, calcula el parámetro Winf y su coeficiente de variación. Por defecto Wfijo=F.}
  \item{Kwi}{Valor de partida de ("start value" o "semilla") del Kwi (parámetro de curvatura) que usa el proceso para iniciar el ajuste del crecimiento en peso.}
  \item{T0wi}{Valor de partida de ("start value" o "semilla") del T0i (punto en el tiempo, en el que el pez tiene un peso cero) que usa el proceso para iniciar el ajuste del crecimiento en peso.}
  \item{n}{Número de réplicas bootstrap. El número recomendado es 1000.}
}
\details{
   Es importante que los valores de partida ("start value" o "semilla") no sean muy diferentes de los que se quiere calcular (se pueden consultar valores de otros años o zonas), si no el proceso
   iterativo por el que se realiza el ajuste podría fallar.

   El ajuste de la ecuación de crecimiento en peso de von Bertalanffy presenta problemas en muchos casos (ver explicación más detallada en Quinn & Deriso (1999)).
   }


\value{
   Dos archivo "txt" con los resultados de los datos originales y de las réplicas bootstrap y su coeficiente de variación. Un archivo contiene los resultados del crecimiento en talla y otro archivo los resultados del crecimiento en peso.


   Dos archivos "pdf" con las gráficas del ajuste con datos originales y bootstrap, las gráficas del análisis de resíduos del ajuste de los datos originales y análisis de las réplicas bootstrap. Un archivo contiene las gráficas del crecimiento en talla y otro archivo las gráficas del crecimiento en peso.


   Si el proceso se realiza separadamente para machos y hembras, los archivos de salida también son independientes.

}

\author{Paz Sampedro}

\references{
  Davison, A.C. & D.V. Hinkley. 1999. Bootstrap methods and their application. Cambridge University Press, 575 pp.

  Efron, B. 1979. Bootstrap methods: another look at the jackknife. Annals of Statistics, 7. 1-26 p.

  Gros, P., Laurec, A. & H. Dupouy. 1987. Estimation des parameters de croissance apport des techniques de reechantillonnage (Bootstrap & Jackniffe) et des statistiques robustes (Médianes et modes).  ICES Comité des statistiques,  CM 1987/D: 10, 25pp.

  Meyer, J.S., Ingersoll, C.G., McDonald, L.L. & M.S. Boyce. 1986. Estimating uncertainty in population growth rates: Jackknife vs Bootstrap techniques. Ecology, 67(5): 1156-1166.

  Moguedet, P. & N. Pérez. Estimation of Megrim (Lepidorhombus whiffiagonis) growth parameters, for males and females, from the ICES division VII: fitting to the V. Bertalanffy model using, resampling techniques, as well as several adjustable central values (mean, median and mode length at age). ICES CM 1988/G: 9. 18 pp.

  Quinn,T.J. & R.B. Deriso. 1999. Quantitative Fish Dynamics. Oxford University Press, 542 pp.
}


\examples{
## Requiere que exista el archivo "bosci2004.csv" en el directorio de trabajo
## Calcula los parámetros de la ecuación de crecimiento en talla y peso de de von Bertalanffy 
## y su coeficiente de variación conjuntamente para ambos sexos

   crecimiento(especie="bosci2004", Li = 40, Ki = 0.27, T0i = -1.9, Wi = 2500, Kwi = 0.4, T0wi = -1)

## Requiere que exista el archivo "rape.csv" en el directorio de trabajo
## Calcula los parámetros K y t0 de la ecuación de crecimiento en talla 
## y su coeficiente de variación; y los parámetros Winf, K y to de la ecuación de crecimiento en peso 
## y su coeficiente de variación
 
   crecimiento(especie="rape", Li= 120, Lfija=T, Ki=0.08, Toi=-2, Wi=10000)
}



