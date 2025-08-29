## FUNCIONES PARA EL PREPROCESADO:

#FUNCION: a traves de una carpeta, unimos todos los individuos de una simulacion en un dataframe.

creacion_df_simulacion <- function(nombre_carpeta){
  df_lista <- list()
  lista_archivos <- list.files(path = paste0("C:/Users/croca/OneDrive/Escritorio/master/TFM_project/TFM_project/TFM/Am66/",nombre_carpeta,"/tracorrTh35/"),full.names = TRUE)
  
  for (i in seq_along(lista_archivos)){
    
    archivo <- read_table(lista_archivos[[i]])
    
    colnames(archivo) <- c('Num_fotogramas', 'X', 'Y')
    
    archivo <- archivo |> 
      mutate(Individuo = sprintf("%03d", i))|>
      mutate(Tiempo =  Num_fotogramas / (25))        #pasamos el num de fotogramas a tiempo (seg))
    
    df_lista[[i]] <- archivo
  }
  
  #unimos los data frames
  data_archivos <- bind_rows(df_lista)|>
    mutate(Individuo = as.factor(Individuo))
  
  #limpiamos los datos: quitar tiempo inicial y final
  data_archivos_limpio <- data_archivos %>% 
    quitar_tiempo_inicial() %>% 
    quitar_tiempo_final() %>% 
    arrange(Individuo,Tiempo) 
  
  
  
  return(data_archivos_limpio)  
}


################################################################################

# FUNCION: quitamos el tiempo inicial en el que aun no ha empezado la evacuacion

quitar_tiempo_inicial <- function(data_archivos){
  
  df_puerta <- marcar_cuadrado_df(data_archivos,longX = 1,longY = 0.5) %>% 
    arrange(Tiempo)
  
  primer_valor <- head(df_puerta$Tiempo, 1)
  individuos_iniciales <- df_puerta %>% 
    filter(Tiempo == primer_valor) %>% 
    filter(cuadrado == 'Si')  %>% 
    pull(Individuo)
  
  for (t in unique(df_puerta$Tiempo)){
    df_tiempo <- df_puerta %>% 
      filter(Tiempo == t) %>% 
      filter(Individuo %in% (individuos_iniciales)) 
    
    if (length(df_tiempo$Individuo) < (length(individuos_iniciales))){
      t_inicio <- t
      break
    }
  }
  data_archivos <- data_archivos %>% 
    filter(Tiempo >= t_inicio) 
  
  return(data_archivos)
}

################################################################################

# FUNCION: quitamos el tiempo final en el que ya queda muy poca gente.

quitar_tiempo_final <- function(data_archivos){
  data_archivos_pers <- data_archivos %>% 
    group_by(Tiempo) %>% 
    summarise(Num_peatones = n())
  data_archivos_pers <- data_archivos_pers %>% 
    mutate(Tiempo_bueno = (ifelse(Num_peatones > max(data_archivos_pers$Num_peatones)*0.2,'Si','No'))) %>% 
    filter(Tiempo_bueno == 'Si') %>% 
    dplyr::select(-c(Tiempo_bueno, Num_peatones))
  
  data_archivos_red <- data_archivos %>% 
    inner_join(data_archivos_pers,'Tiempo') 
  
  return(data_archivos_red)
}


#############################################################################

#FUNCION: crear df de todas las simulaciones a traves de una lista de archivos y marcarlos con el id
creacion_df_total_simulaciones <- function(lista_archivos){
  
  # Generar df de cada simulacion y guardarlos con el id
  df <- list()
  
  for (i in seq_along(lista_archivos)){
    
    archivo <- lista_archivos[[i]]
    
    df[[i]] <- creacion_df_simulacion(nombre_carpeta = archivo) %>% 
      mutate(id = substr(archivo,nchar(archivo) - 5, nchar(archivo)))
  }
  data_archivos_total <- bind_rows(df)
  
  return(data_archivos_total)
}

################################################################################

#FUNCION: nos quedamos con el espacio de estudio a un cuadrado

individuos_cuadrado_df <- function(data_archivos, longX, longY) {
  df_cuadrado <- data_archivos |> 
    mutate(cuadrado = ifelse(abs(X) <= longX & Y <= longY, 'Si', 'No')) |> 
    filter(cuadrado == 'Si') |> 
    dplyr::select(-cuadrado)  # Asegurar que select() es de dplyr
  
  return(df_cuadrado)  
}


################################################################################
################################################################################


# FUNCIONES CROCKER PLOT

cortar_data <- function(data,n_breaks){
  
  #ordenamos
  data <- data %>% arrange(Num_fotogramas)
  
  # Obtener valores únicos de Num_fotogramas, 
  valores_unicos <- unique(data$Num_fotogramas)
  
  # Dividir esos valores únicos en n_breaks grupos
  grupos_unicos <- cut(seq_along(valores_unicos), breaks = n_breaks, labels = FALSE)
  
  # Crear tabla de asignación
  asignacion <- data.frame(
    Num_fotogramas = valores_unicos,
    grupo = grupos_unicos
  )
  
  # Unir con el data original
  data_dividida <- data %>%
    left_join(asignacion, by = "Num_fotogramas")
  
  # Dividir en lista
  lista_df <- data_dividida %>%
    group_split(grupo)
  
  return(lista_df)
}

################################################################################

# FUNCION: dividimos el data en distintos partes segun el num de fotogramas que queremos analizar (n_t)
#corremos uno porque luego haremos saltos de 5 entonces cada tiempo sera distinto
datas_consecutivos <- function(data, n_t) {
  
  lista_df <- list()
  data <- data %>% arrange(Num_fotogramas)
  valores_unicos <- unique(data$Num_fotogramas)
  
  for (i in 1:(length(valores_unicos) - n_t + 1)) {
    fotogramas_seleccionados <- valores_unicos[i:(i + n_t - 1)]
    
    data_red <- data %>%
      filter(Num_fotogramas %in% fotogramas_seleccionados)
    
    lista_df[[i]] <- data_red
  }
  
  return(lista_df)
}

################################################################################

#FUNCION1
#para cada instante de tiempo, calculo los diagramas de persistencia y los almaceno en homologydata

getIntervals <- function(mydata, L, mymaxdimension, mymaxscale, thistime) {
  # Normaliza las posiciones x e y
  mydata$X <- mydata$X / L
  mydata$Y <- mydata$Y / L
  
  # Calcula la matriz de distancias entre puntos usando las posiciones x e y
  dist_matrix <- as.matrix(dist(cbind(mydata$X, mydata$Y)))
  
  # Calcula el complejo Vietoris-Rips para las distancias
  results <- ripsDiag(X = dist_matrix, dist = "arbitrary", maxdimension = mymaxdimension, 
                      maxscale = mymaxscale, library = "GUDHI", printProgress = FALSE)
  
  # Obtiene el diagrama de persistencia
  results <- results$diagram
  
  # Devuelve los intervalos con la información temporal
  return(data.frame(t = rep(thistime, nrow(results)), dimension = results[, 1], birth = results[, 2], death = results[, 3]))
}

################################################################################

#FUNCION2: crea un df con el numero de betti para cada epsilon y tiempo,
#le pasas el df de la homologia perssistente y el numero de epsilons a generar
turnIntervalsIntoGrid <- function(homologydata, numEpsilons,mymaxscale) {
  dimensions <- unique(homologydata$dimension)
  times <- unique(homologydata$t)
  
  # Secuencia de epsilon
  epsilons <- seq(from = 0, to = mymaxscale, length.out = numEpsilons)
  
  # Resultado final
  bettiData <- vector("list", length(dimensions) * length(times))
  idx <- 1
  
  for (thisDimension in dimensions) {
    for (thisTime in times) {
      
      # Filtrar por dimensión y tiempo
      snapshot <- subset(homologydata, dimension == thisDimension & t == thisTime)
      if (nrow(snapshot) == 0) next
      
      births <- snapshot$birth
      deaths <- snapshot$death
      
      # Vector donde contaremos los betti numbers por epsilon
      counts <- sapply(epsilons, function(eps) {
        sum(births <= eps & deaths >= eps)
      })
      
      # Agregar datos al resultado
      bettiData[[idx]] <- data.frame(
        t = thisTime,
        dimension = thisDimension,
        epsilon = epsilons,
        betticount = counts
      )
      idx <- idx + 1
    }
  }
  
  # Combinar todos los resultados
  bettiData <- do.call(rbind, bettiData[1:(idx - 1)])
  return(bettiData)
}

################################################################################

#FUNCION3: crea el grafico: df de bettiData, la dim que quieres visualizar, el espilon maximo (eje y)
crockerplot_grafic <- function(bettiData, dim, maxYlim) {
  # genera gráficos de nivel para la dimensión dada
  
  df <- subset(bettiData, dimension == dim)
  
  # Calcular los niveles de contorno de forma automática
  betti_range <- range(df$betticount, na.rm = TRUE)
  
  if (dim == 0){
    breaks <- c(0:60)
  }else{
    breaks <- c(0:20)
  }
  
  #calcula curvas de nivel a partir de la variable betticount
  p <- ggplot(df, aes(x = t, y = epsilon, z = betticount)) +
    geom_contour(breaks = breaks, aes(colour = ..level..)) +    #...level... representa los valores de los contornos, 
    theme_bw() +
    scale_colour_gradientn(
      limits = range(breaks),
      colours = rainbow(length(breaks), end = 0.83),
      guide = "legend"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    ylim(0, maxYlim) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    labs(x = "Simulation time t", y = "Proximity Parameter ε") + 
    ggtitle(paste("Crocker plot for dimension", dim))
  
  return(p)
}

################################################################################
#FUNCION: crea un gráfico que muestra cómo cambia el número de agujeros topológicos en el tiempo y según la proximidad entre puntos.

matrizBetti <- function(bettiData, dim, maxYlim) {
  # Filtrar por dimensión
  df <- subset(bettiData, dimension == dim)
  
  p <- ggplot(df, aes(x = t, y = epsilon, fill = betticount)) +
    geom_tile() +
    #geom_text(aes(label = betticount), color = "black", size = 3) + 
    theme_bw() +
    scale_fill_gradientn(
      colours = rainbow(6),
      limits = c(0, max(df$betticount, na.rm = TRUE)),
      name = "Betti count"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    ylim(0, maxYlim) +
    labs(
      x = "Simulation time t",
      y = "Proximity Parameter ε",
      title = paste("Matriz de BettiCounts for Dimension", dim)
    )
  
  return(p)
}
