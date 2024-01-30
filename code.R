#########################
# Estudio de estabilidad de puntos de equilibrio de
# Ecuaciones en diferencias autónomas
# Autor: Diego Mestanza
# Asignatura: Temas de Matemáticas Avanzadas
#########################

# 1. Paquetes necesarios ----

#install.packages("Derive")

library(polynom)
#library(Derive)
library(tidyverse)

library(showtext)
#library(sysfonts)
library(randomcoloR)

showtext_auto()
font_add_google("Montserrat", "Montserrat")

my_theme <- theme_bw() +
  theme(
    panel.grid.major = element_line(colour = "#f0f0f0"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(),
    text = element_text(size = 12, family = "Montserrat"),
    plot.title = element_text(face="bold"),
    plot.subtitle = element_text(size = 14),
    plot.caption = element_text(hjust = 0, face = "italic")
  )

# 2. Función que permite generar el diagrama de fases, las trayectorias y evaluar
# analíticamente una ecuación en diferencias, a partir de condiciones iniciales ----

DE <- function(coef_f,
               coef_f2,
               coef_g,
               two = FALSE,
               min = -6,
               max = 6,
               step = 0.2,
               initials){
  
  if (two) {
    expr_f   <- polynomial(coef = coef_f)
    expr_f2  <- polynomial(coef = coef_f2)
    expr_g   <- polynomial(coef = coef_g)
    
    deriv_f  <- deriv(expr_f, name = c("x"))
    deriv_f2 <- deriv(expr_f2, name = c("x"))
    deriv    <- (deriv_f*expr_f2 - deriv_f2*expr_f)
    deriv2   <- expr_f2^2
    
    roots_g  <- solve(expr_g)
    
    for (root in roots_g) {
      evaluate <- predict(deriv, root)/predict(deriv2, root)
      if (abs(evaluate) < 1) {
        print(paste0(round(root, digits = 2), " es un punto LAE"))
      }
      else {
        print(paste0(round(root, digits = 2), " es un punto inestable"))
      }
    }
    
    x <- seq(min, max, step)
    evaluate <- predict(expr_f, x)/predict(expr_f2, x)
    
  } else {
    expr_f  <- polynomial(coef = coef_f)
    expr_g  <- polynomial(coef = coef_g)
    
    deriv   <- deriv(expr_f, name = c("x"))
    
    roots_g <- solve(expr_g)
    
    for (root in roots_g) {
      evaluate <- predict(deriv, root)
      if (abs(evaluate) < 1) {
        print(paste0(round(root, digits = 2), " es un punto LAE"))
      } else {
        print(paste0(round(root, digits = 2), " es un punto inestable"))
      }
    }
    
    x <- seq(min, max, step)
    evaluate <- predict(expr_f, x)
    
  }
  
  y <- x
  
  data <- data.frame(x, y, evaluate)

  for (initial in initials) {
    
    data <- data %>% 
      mutate(cx = 0, cy = 0, xt = 0, time = seq(0, length(data$x)-1, 1))
    
    data$cx[1:2] <- initial
    data$xt[1] <- initial
    
    for (p in seq(2, length(data$x), 2)) {
      
      if (two) {
        data$cy[p]   <- predict(expr_f, data$cx[p])/predict(expr_f2, data$cx[p])
      } else {
        data$cy[p]   <- predict(expr_f, data$cx[p])
      }
      
      data$cx[p+1] <- data$cy[p]
      data$cy[p+1] <- data$cx[p+1]
      if (p != length(data$x)-1) data$cx[p+2] <- data$cy[p+1]
      
    }   
    
    for (p in seq(2, length(data$x), 1)) {

      if (two) {
        data$xt[p]   <- predict(expr_f, data$xt[p-1])/predict(expr_f2, data$xt[p-1])
      } else {
        data$xt[p]   <- predict(expr_f, data$xt[p-1])
      }      
    }
    
    if (initial == initials[1]) {
      
      plot_fases <- ggplot(data) +
        geom_line(aes(x = x, y = evaluate), 
                  colour = "blue", 
                  linewidth = 1.4) +
        geom_line(aes(x = x, y = y), 
                  colour = "black", 
                  linewidth = 1.4) +
        geom_point(aes(x = cx, y = cy),
                   alpha = 0,
                   na.rm = TRUE) +
        geom_path(aes(x = cx, y = cy), 
                  colour = randomcoloR::randomColor(),
                  linewidth = 1,
                  na.rm = TRUE) +
        labs(y = "", x = "", title = "Diagrama de Fases") + 
        coord_cartesian(xlim = c(min, max), ylim = c(min, max)) +
        my_theme
      
      plot_trajectory <- ggplot(data) +
        geom_point(aes(x = time, y = xt),
                   alpha = 0,
                   na.rm = TRUE) +
        geom_path(aes(x = time, y = xt), 
                  colour = randomcoloR::randomColor(),
                  linewidth = 1,
                  na.rm = TRUE) +
        labs(y = "", x = "", title = "Trayectorias") +
        my_theme
      
    } else {
      
      plot_fases <- plot_fases +
        geom_point(data = data, aes(x = cx, y = cy), 
                   alpha = 0, 
                   na.rm = TRUE) +
        geom_path(data = data, aes(x = cx, y = cy), 
                  colour = randomcoloR::randomColor(),
                  linewidth = 1, 
                  na.rm = TRUE)
      
      plot_trajectory <- plot_trajectory +
        geom_point(data = data, aes(x = time, y = xt), 
                   alpha = 0, 
                   na.rm = TRUE) +
        geom_path(data = data, aes(x = time, y = xt), 
                  colour = randomcoloR::randomColor(),
                  linewidth = 1, 
                  na.rm = TRUE)
      
    }
  }
  print(plot_fases)
  print(plot_trajectory)
}



# 3. Evaluar la función generada en 2. con la ecuación en diferencias y 
# las condiciones iniciales propuestas en la tarea de Temas de Matemáticas Avanzadas ----

DE(coef_f = c(-2, 2, 1, 1), 
   coef_f2 = c(1, 0, 1), 
   coef_g = c(2, -1, -1), 
   two = TRUE, 
   min = -6,
   max = 6,
   step = 0.2,
   initials = c(-6, 0.9, 2))
