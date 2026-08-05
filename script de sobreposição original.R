#Para script das análises de nicho serão necessárias:
#  1- Dados de ocorrência filtrados das espécies.
#  2- Dados de Plano de fundo (Background).
#  3- Variáveis ambientais recortadas.

# Carregamento dos pacotes----
# Aqui ta baixando um pacote chamado pacman, que carrega e baixa arquivos de uma vez e mais rápido, caso não tenha rode a linha abaixo sem o # 
#if(!require(pacman)) install.packages('pacman', repos = c('https://trinker.r-universe.dev', 'https://cloud.r-project.org'))
pacman::p_load(sp,ecospat,here,tidyverse,sf,randomcoloR,stars,readxl,colorspace,geodata,ade4,terra,reshape2,knitr)
getwd()
dir()
# Carregamento dos pontos  de  Ocorrência----
#Crie em seu diretório uma pasta para seus Dados, que terá seus pontos de ocorrência que serão filtrados 
#para isso clique aqui ao lado na aba FILES aperte em NEW FOLDER, nesse script minha pasta se chama "Dados" mas mude como desejar
#caso queira usar o R rode abaixo sem o #
## Criar diretório----
dir.create(here::here("Dados"))  #Pasta dos Dados Poõem aqui as ocorrẽncias (Tem que ta no formato "Especie\Longitude\Latitude)
dir.create(here::here("Imagens")) #Pasta das Imagens
dir.create(here::here("Tabelas"))  #Pasta das Tabelas
dir.create(here::here("Preditores"))  #Pasta dos Preditores recortados

#Carregar os pontos----
pontos_occ <- read_excel(here("Dados",
                              "Alouatta_Completo.xlsx"))
pontos.limpos<-na.omit(pontos_occ)

Especie <- unique(pontos.limpos$Especie) %>%                                     #Checagem das ocorências
  print()


pontos_occ.list <- pontos_occpontos_occ.list <- list()                                                #Fazer uma lista com os pontos de ocorrencia
for (i in seq_along(Especie)) {
  pontos_occ.list[[i]] <- as.data.frame(pontos.limpos[pontos.limpos$Especie == Especie[i], 2:3])
}
names(pontos_occ.list) <- Especie

#Definir a quantidade de grupos analisados----
n.groups <- 2
g.names <- Especie
g.codenames <- g.names
g.colors <- distinctColorPalette(5) #cores aleatórias igual a quantidade de espécies - Opção automática
g.colors <- c("#C906CF","#16941E") 

#Checagem dos pontos de ocorrência no mapa


pontos_occ.list2 <- do.call(rbind, pontos_occ.list) %>%
  tibble::rownames_to_column(var = "Especie") %>%
  mutate(Especie = gsub("\\..*", "", Especie),
         Especie = factor(Especie, g.names))

mapa <- sf::st_read(dsn ="C:/Users/WORKSTATION/Igor/Neotropico/Am_Sul_Completa.shp" ,
                    quiet = T)# Parte que escolhe o Mapa, substituir pelas áreas de endemismo
sf::sf_use_s2(FALSE)

limits <- range(pontos_occ.list2$longitude)
lat_limits <- range(pontos_occ.list2$latitude)

g <- ggplot() +                                                    #Verificação dos pontos no mapa
  geom_sf(data = mapa) +  
  theme_minimal() +
  geom_point(pontos_occ.list2, mapping = aes(longitude, latitude,
                                             col = Especie),
             size = 3) + 
  scale_color_manual(values = g.colors) +  
  xlab("") +
  ylab("") +
  coord_sf(xlim = limits, ylim = lat_limits) +
  
  theme(
    strip.text = element_text(face = "italic"),
    legend.text = element_text(face = "italic"),
    legend.position = "none"
  )
g
loc <- here("Imagens", "Mapa ocorrências Alouatta.tiff") # para modificar o nome da pasta, 
ggsave(loc, g)

# Filtrar para os países que estão localizados

map <- mapa[mapa$nome %in% c("Brasil"), ]  
pnts_sf <- st_as_sf(x = pontos_occ.list2,
                    coords = c('longitude', 'latitude'),
                    crs = st_crs(map))

pnts <- pnts_sf %>%
  mutate(intersection = as.integer(st_intersects(geometry, map)))

remove <- is.na(pnts$intersection)
pontos_occ.list3 <- pontos_occ.list2[!remove, ]

#Escolher a forma de Área de calibração----
#Nessa versão estamos utilizando um Mínimo Polígono Convexo (MPC) que é o original do Bruno Vilela, 
#mas vamos testar depois de maneiras diferentes. Foi usado 2º com base na dispersão da espécie

buffer.size <- 1                             #Lembrando que 1° -> 111 km
mcp <- function(xy) {
  hull <- xy %>% 
    st_as_sf(coords = c("longitude", "latitude")) %>% 
    st_union() %>% 
    st_convex_hull()
  return(hull)
}

#Carregamentodos dos Dados ambientais----
variavel <- list.files (path = "Preditores/Current", pattern = "\\.tif$", full.names = T) %>% #Caminho para a pasta das variáveis recortadas
  terra::rast()

names(variavel)                                                    #Checagem visual das  variáveis
plot(variavel)

#Preparação dos Dados----
#Agora que temos as 3 partes vamos juntar para analisar! Aqui embaixo vamos usar os pontos para gerar um Mínimo Polígono Convexo 
#com o buffer da área de calibração. Após isto veremos os valores dos Dados ambientais por grupo dos pontos de ocorrência + área de calibração. 
#Por fim vamos visualizar em um plot

# Objetos nulos
g.assign <- pontos_occ.list3$Especie
xy.mcp <- list()
back.env <- list()
spec.env <- list()
row.sp <- list()

united <- st_union(st_make_valid(map))

# Loop
for (i in 1:n.groups) {
  # Salvar as linhas por nº de espécie
  g.limit <- g.assign == Especie[i]
  row.sp[[i]] <- which(g.limit)
  
  # Polígono de fundo 
  mcp.occ <- mcp(pontos_occ.list3[g.limit, -1])
  xy.mcp.i <- st_buffer(mcp.occ, dist = buffer.size)
  st_crs(xy.mcp.i) <- st_crs(united)
  xy.mcp[[i]] <-  st_as_sf(st_intersection(xy.mcp.i, united))
  xy.mcp[[i]]$Especie <- Especie[i]
  # Ambiente de calibração
  extract_temp <- terra::extract(variavel, xy.mcp[[i]])[, -1]
  back.env[[i]] <- na.exclude(extract_temp)
  # Ambiente das espécies
  spec.env[[i]] <- na.exclude(terra::extract(variavel, 
                                             pontos_occ.list3[g.limit, -1])[, -1])
}

#Buffer nos mapas:
xy.mcp.is <- do.call(rbind, xy.mcp)
xy.mcp.is$Especie <- factor(xy.mcp.is$Especie, xy.mcp.is$Especie)


#Ver os mapas  com os buffers
g <- ggplot() +
  geom_sf(data = map) + #Atenção, nessa parte filter(world, ???) é aonde que o plot do R 
  theme_minimal() +                                                         #vai dar o zoom no mapa, seja ele em um continente (region_wb) ou país (iso_a3)
  geom_sf(data = xy.mcp.is, aes(col = Especie)) +
  geom_point(pontos_occ.list3, mapping = aes(longitude, latitude,
                                             col = Especie),
             size = .3) +
  scale_color_manual(values = g.colors) +
  xlab("") +
  ylab("") +
  theme(
    strip.text = element_text(face = "italic"),
    legend.text = element_text(face = "italic"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1))
g
loc <- here("Imagens", "Mapa buffer Alouatta.tiff")
ggsave(loc, g)

#Salvar os pontos de ocorrência numa tabela.csv
write.csv(pontos_occ.list3, 
          file = here("Dados", "ocorrências Alouatta.csv"),
          row.names = FALSE)


#Organizar as tabelas que serão utilizadas

# Pontos de ocorrência por grupo
g.occ.points <- pontos_occ.list3
colnames(g.occ.points)[1] <- "Groups"
# Valores ambientais da área de calibração 
all.back.env <- do.call(rbind.data.frame, back.env)
# Valores ambientais para cada ponto de ocorrência
all.spec.env <- do.call(rbind.data.frame, spec.env)
# Valores ambientais de todas as  espécies juntas
data.env <- rbind(all.spec.env, all.back.env)

#Checagem de números de pontos de ocorrência por grupo
table(g.occ.points[, 1])

#Comparação de nicho----
#As análises de nicho seguem estes passo a passo
# 1 Análise das variáveis ambientais por meio de uma PCA 
# 2 Medição da Sobreposição de Nicho por índice de Schoener "D"
# 3 Ver as dinâmicas de nicho: Não ocupado (Unfilling), Estável (Stability), Em expansão (Expansion)

# 1º Cálculo da PCA----
# Matriz de peso
w <- c(rep(0, nrow(all.spec.env)), rep(1, nrow(all.back.env)))
# PCA das variáveis ambientais
pca.cal <- dudi.pca(data.env, row.w = w, center = TRUE, 
                    scale = TRUE, scannf = FALSE, nf = 2)

#Após termos os resultados da PCA, vamos utilizar o valor do primeiro e segundo vetor (eigenvector) 
#dela para os valores de calibração e ocorrência de cada grupo

# Linha da planilha que correspondem aos valores de sp1
adtion <- cumsum(c(0, sapply(back.env, nrow)))
begnd <- nrow(all.spec.env)
# Criação de listas vazias para salvar os resultados
scores.back <- list()
scores.spec <- list()

# Adicionar os valores 
for (i in 1:n.groups) {
  scores.spec[[i]] <- pca.cal$li[row.sp[[i]], ]
  pos <- (begnd[1] + adtion[i] + 1) : (begnd[1] + adtion[i + 1])
  scores.back[[i]] <- pca.cal$li[pos, ]  
}

total.scores.back <- do.call(rbind.data.frame, scores.back)

#Espaço ambiental
#Aqui vamos gerar o espaço ambiental com base no valor da PCA que acabamos de calcular utilizando 
#a área de calibração e os registros de ocorrência.

#Definição da resolução do grid do Espaço Ambiental
R <- 100

#Em seguida, modelamos a densidade das espécies no grid ambiental, considerando a quantidade de espécies por var.ambiental da área de calibração
z <- list()

for (i in 1:n.groups) {
  z[[i]] <- ecospat.grid.clim.dyn(total.scores.back,
                                  scores.back[[i]],
                                  scores.spec[[i]],
                                  R = R)
}

# 2º Medição da Sobreposição de Nicho----
#O valor da Sobreposição de Nicho é dado pelo D de Schoener, que varia de 0 a 1, 
#quanto mais a sobreposição mais próximo de 1 é o resultado.

#O D de Schoener é calculado por um teste de similaridade, aqui abaixo escolhemos quantas repetições vamos realizar.
#Este numero de repetições será usado nos valores de análises de dinâmicas de nicho
rep <- 100

#Após o numero de repetições escolhidos, podemos gerar os valores. 
# Matrizes Vazias
D <- matrix(nrow = n.groups, ncol = n.groups)
rownames(D) <- colnames(D) <- g.codenames
unfilling <- stability <- expansion <- equi <- sim <- D


for (i in 2:n.groups) {

  for (j in 1:(i - 1)) {
    cat(sprintf("Calculando grupo %d vs grupo %d...\n", j, i))
    x1 <- z[[j]]
    x2 <- z[[i]]
    
    # Resultados do D Schoener e Dinâmicas de Nicho
    
    # Sobreposição de Nicho
    D[i, j] <- ecospat.niche.overlap(x1, x2, cor = TRUE)$D
    
  
    # Similaridade de Nicho (valor de p)
    sim[i, j] <- ecospat.niche.similarity.test(x1, x2, rep)$p.D
    sim[j, i] <- ecospat.niche.similarity.test(x2, x1, rep)$p.D
    
    # Equivalência de Nicho (valor de p)
    equi[i, j] <- ecospat.niche.equivalency.test(x1, x2, rep)$p.D
    equi[j, i] <- ecospat.niche.equivalency.test(x2, x1, rep)$p.D
    
    
    
    # Nichos Em expansão (expansion), Estáveis (stability) e Não ocupados (unfilling).  
    index1 <- ecospat.niche.dyn.index(x1, x2)$dynamic.index.w
    index2 <- ecospat.niche.dyn.index(x2, x1)$dynamic.index.w
    expansion[i, j] <- index1[1]
    stability[i, j] <- index1[2]
    unfilling[i, j] <- index1[3]
    expansion[j, i] <- index2[1]
    stability[j, i] <- index2[2]
    unfilling[j, i] <- index2[3]
  }
}

## Resultados numéricos

#Valor de D de Schoener:
kable(D, digits = 3, format = "markdown")

round(D, digits = 3)

write.table(D,
            file = here::here("Tabelas", "D de Schoener Alouatta.txt"),
            sep = " ")



D.teste <- D 
D.teste <- round(D.teste, digits = 3)
D.teste <- melt(D.teste)


D_heatmap <- ggplot(D.teste, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = value), size = 4) +
  scale_fill_gradient2(low  = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 0, na.value = "grey92",
                       name = "Valores de Sobreposição") +
  labs(title    = "D de Schoener",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
D_heatmap

loc <- here("Imagens", "D de Schoener Alouatta.tiff") # para modificar o nome da pasta, 
ggsave(loc, D_heatmap) 





#Modelo nulo de similaridade de nicho (valor de p):

kable(sim, digits = 3, format = "markdown")

round(sim, digits = 3)


sim.teste <- sim 
sim.teste <- round(sim.teste, digits = 3)
sim.teste <- melt(sim.teste)


Sim_heatmap <- ggplot(sim.teste, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = value), size = 4) +
  scale_fill_gradient2(low  = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 0, na.value = "grey92",
                       name = "Valores de p") +
  labs(title    = "Valor de P (teste de similaridade)",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
Sim_heatmap

loc <- here("Imagens", "Valor de p (teste sim) Alouatta.tiff") # para modificar o nome da pasta, 
ggsave(loc, Sim_heatmap) 


cbind("Espécies",sim)
write.table(sim,
            file = here::here("Tabelas", "sim_Alouatta.txt"),
            sep = " ")




#Equivalência de nicho:

kable(equi, digits = 3, format = "markdown")

round(equi, digits = 3)


equi.teste <- equi 
equi.teste <- round(equi.teste, digits = 3)
equi.teste <- melt(equi.teste)


Equi_heatmap <- ggplot(equi.teste, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = value), size = 4) +
  scale_fill_gradient2(low  = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 0, na.value = "grey92",
                       name = "Equivalência de nicho") +
  labs(title    = "Valor de P Equivalência",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
Equi_heatmap

loc <- here("Imagens", "Equivalência_Alouatta.tiff") # para modificar o nome da pasta, 
ggsave(loc, Equi_heatmap) 


cbind("Espécies",equi)
write.table(equi,
            file = here::here("Tabelas", "equi_Alouatta.txt"),
            sep = " ")


#Nichos Não Ocupados:

kable(unfilling, digits = 3,  format = "markdown")

round(unfilling, digits = 3)


unfilling.teste <- unfilling 
unfilling.teste <- round(unfilling.teste, digits = 3)
unfilling.teste <- melt(unfilling.teste)


unf_heatmap <- ggplot(unfilling.teste, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = value), size = 4) +
  scale_fill_gradient2(low  = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 0, na.value = "grey92",
                       name = "Valores de Nichos") +
  labs(title    = "Nichos Perdidos",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
unf_heatmap

loc <- here("Imagens", "Unfilling_Alouatta.tiff") # para modificar o nome da pasta, 
ggsave(loc, unf_heatmap) 



cbind("Espécies",unfilling)
write.table(unfilling,
            file = here::here("Tabelas", "unfilling_Alouatta.txt"),
            sep = " ")


#Nichos em Expansão:

kable(expansion, digits = 3,  format = "markdown")

round(expansion, digits = 3)

expansion.teste <- expansion 
expansion.teste <- round(expansion.teste, digits = 3)
expansion.teste <- melt(expansion.teste)


Exp_heatmap <- ggplot(expansion.teste, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = value), size = 4) +
  scale_fill_gradient2(low  = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 0, na.value = "grey92",
                       name = "Valores de Expansão") +
  labs(title    = "Nicho em expansão",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
Exp_heatmap

loc <- here("Imagens", "Expansion_Alouatta.tiff") # para modificar o nome da pasta, 
ggsave(loc, Exp_heatmap) 


cbind("Espécies",expansion)
write.table(expansion,
            file = here::here("Tabelas", "expansion_Alouatta.txt"),
            sep = " ")


#Nichos Estáveis:

kable(stability, digits = 3,  format = "markdown")

round(stability, digits = 3)


stability.teste <- stability 
stability.teste <- round(stability.teste, digits = 3)
stability.teste <- melt(stability.teste)


Est_heatmap <- ggplot(stability.teste, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = value), size = 4) +
  scale_fill_gradient2(low  = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 0, na.value = "grey92",
                       name = "Valores de Nicho") +
  labs(title    = "Nicho estáveis",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
Est_heatmap

loc <- here("Imagens", "Stability_Alouatta.tiff") # para modificar o nome da pasta, 
ggsave(loc, Est_heatmap) 


cbind("Espécies",stability)
write.table(stability,
            file = here::here("Tabelas", "stability_Alouatta.txt"),
            sep = " ")






### Imagens (plot) dos diversos nichos
#Atenção nessa parte não consegui ainda
#plot base:
for (i in 2:n.groups) {
  
  for (j in 1:(i - 1)) {
    
    arquivo <- here::here( "Imagens", paste0("Nicho_", Especie[[1]], "_x_", Especie[[2]], ".png"))
    
    png(
      filename = arquivo,
      width = 2000,
      height = 2000,
      res = 300
    )
    
    
    col1 <- colorRampPalette(c(desaturate(g.colors[i]), g.colors[i]))(5)
    col2 <- colorRampPalette(c(desaturate(g.colors[j]), g.colors[j]))(5)
    col_int <- colorRampPalette(c(desaturate('#e7298a'), '#e7298a'))(5)
    Nich.din<- ecospat.plot.niche.dyn(z[[i]], z[[j]], 0, 
                                      colZ1 = g.colors[i],
                                      colZ2 = g.colors[j],
                                      colinter = adjustcolor(col_int, .4),
                                      colz1 = adjustcolor(col1, .4),
                                      colz2 = adjustcolor(col2, .4),
                                      name.axis1 = "PC1",
                                      name.axis2 = "PC2") + mtext(paste(Especie[[1]],"x", Especie[[2]]), font = 3)
    dev.off()
     }
}



arquivo <- here::here( "Imagens", paste0("Nicho_", Especie[[1]], ".png"))

png(
  filename = arquivo,
  width = 2000,
  height = 2000,
  res = 300
)

Nicho.esp <- ecospat.plot.niche (z[[1]], 
                                 title = "Nicho Alouatta Amazônia",
                                 name.axis1 = "PC1",
                                 name.axis2 = "PC2",
                                 cor=FALSE)
dev.off()
ggsave(loc, Nicho.esp) 





arquivo <- here::here( "Imagens", paste0("Nicho_", Especie[j], "_x_", Especie[[2]], ".png"))
png(
  filename = arquivo,
  width = 2000,
  height = 2000,
  res = 300
)

Nicho.esp.2 <- ecospat.plot.niche (z[[2]], 
                                   title = "Nicho Alouatta Mata Atlântica",
                                   name.axis1 = "PC1",
                                   name.axis2 = "PC2",
                                   cor=FALSE)
dev.off()
ggsave(loc, Nicho.esp.2) 


Nicho.esp <- ecospat.plot.niche (z[[1]], 
                                 title = "Área Alouatta Amazônia ",
                                 name.axis1 = "PC1",
                                 name.axis2 = "PC2",
                                 cor=FALSE)

dev.off()
ggsave(loc, Nicho.esp) 




#Imagem da contribuição de cada variável para o eixo da PCA. 
#Verifique o código da variável em <http://www.worldclim.org/bioclim> ou em sua planilha.

loadings <- cbind(cor(data.env, pca.cal$tab[,1]), cor(data.env, pca.cal$tab[,2]))
colnames(loadings) <- c("axis1", "axis2")
loadings <- loadings[c(1, 12:19, 2:11), ]

barplot(loadings[,1], las=2, main="PC1")

barplot(loadings[,2], las=2, main="PC2")

#Visualização em setas diretamente no espaço ambiental.
contrib <- pca.cal$co
eigen <- pca.cal$eig
nomes <- numeric(20)
for(i in 1:20){
  nomes[i] <- paste('bio',i, sep="")
}
s.corcircle(contrib[, 1:2] / max(abs(contrib[, 1:2])), 
            grid = F,  label = nomes, clabel = 1.2)
text(0, -1.1, paste("PC1 (", round(eigen[1]/sum(eigen)*100,2),"%)",
                    sep = ""))
text(1.1, 0, paste("PC2 (", round(eigen[2]/sum(eigen)*100,2),"%)",
                   sep = ""), srt = 90)




#####testes pos analise

str(index1)
index1


getAnywhere(ecospat.niche.dyn.index)

index1 <- ecospat.niche.dyn.index(x1, x2)

index1$category_quantity
save.image(file =  "Alouatta.RData" )
