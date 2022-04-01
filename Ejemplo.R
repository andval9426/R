
naive_results <- import_results(file="Datos/pubmed-FundacionV-set.nbib")

# Cargar librerias ----
library(dplyr); library(ggplot2); library(ggraph); library(igraph)
library(readr); library(devtools); library(easyPubMed); library(RISmed)

# Instalar paquete ----
#install_github("elizagrames/litsearchr", ref="main")
library(litsearchr)
packageVersion("litsearchr")


# Cargar resultados pubmed ----
naive_results <- import_results(file="Datos/pubmed-medication-set.nbib")
naive_results2 <- import_results(file="Datos/scopus.bib")

# 1. Cargar resultados a partir de APIs PUBMED ----
my_query <-  "(medication) AND (CBT) AND (phobia)"
my_entrez_id <- get_pubmed_ids(my_query)
my_abstracts_xml <- fetch_pubmed_data(pubmed_id_list = my_entrez_id)
my_titles <- custom_grep(my_abstracts_xml, "ArticleTitle", "char")
TTM <- nchar(my_titles) > 75
my_titles[TTM] <- paste(substr(my_titles[TTM], 1, 70), "...", sep = "")
head(my_titles,10)

# Explorando la busqueda realizada en PUBMED ----
nrow(naive_results) # numero de filas 
colnames(naive_results) # nombres de las columnas
naive_results[5, "keywords"] # palabras claves fila 5
sum(is.na(naive_results[, "keywords"])) # numero de paper sin palabras claves
extract_terms(keywords=naive_results[, "keywords"],
              method="tagged") # extraer palabras claves


# palabras claves ----
keywords <- extract_terms(keywords=naive_results[, "keywords"], 
                          method="tagged", min_n=1)

# palabras vacias ----
clinpsy_stopwords <- read.table("Datos/palabras_vacias.txt", header = TRUE)
clinpsy_stopwords <- clinpsy_stopwords$palabras_vacias

# Union palabras vacias en ingles y las anteriores ----
all_stopwords <- c(get_stopwords("English"), clinpsy_stopwords)

# Extraer solo los términos de busqueda relevantes ----
title_terms <- extract_terms(
  text=naive_results[, "title"],
  method="fakerake",
  min_freq=3, min_n=2,
  stopwords=all_stopwords
)

title_terms

# Terminemos sumando los términos de búsqueda que obtuvimos de los títulos de 
# los artículos y los que obtuvimos de las palabras clave anteriormente, eliminando 
# los duplicados.
terms <- unique(c(keywords, title_terms))


# Analisis de RED ----
docs <- paste(naive_results[, "title"], naive_results[, "abstract"])
docs[1]

# Creacion de matriz de documentos ----
dfm <- create_dfm(elements=docs, features=terms)
g <- create_network(dfm, min_studies=3)

# Visualizacion de la red ----
ggraph(g, layout="stress") +
  coord_fixed() +
  expand_limits(x=c(-3, 3)) +
  geom_edge_link(aes(alpha=weight)) +
  geom_node_point(shape="circle filled", fill="white") +
  geom_node_text(aes(label=name), hjust="outward", check_overlap=TRUE) +
  guides(edge_alpha=FALSE)


# PODA ----
## La 'fuerza' de cada término en la red es el número de otros términos con los que aparece junto
strengths <- strength(g)
data.frame(term=names(strengths), strength=strengths, row.names=NULL) %>%
  mutate(rank=rank(strength, ties.method="min")) %>%
  arrange(strength) ->
  term_strengths

term_strengths

## Grafico de la "fuerza": para eliminar los términos de búsqueda menos vinculados a los demás.
cutoff_fig <- ggplot(term_strengths, 
                     aes(x=rank, y=strength, label=term)) +
  geom_line() +
  geom_point() +
  geom_text(data=filter(term_strengths, rank>5), hjust="right", nudge_y=20, check_overlap=TRUE)

cutoff_fig

## Corte
cutoff_cum <- find_cutoff(g, method="cumulative", percent=0.8)

cutoff_cum

## Representacion basandose en el corte
cutoff_fig +
  geom_hline(yintercept=cutoff_cum, linetype="dashed")

## Una vez que hemos encontrado un valor de corte, 
## la reduce_graph()función lo aplica y elimina los términos con poca fuerza. 
## Los argumentos son la red original y el corte. 
## La get_keywords()función a continuación, obtiene los términos restantes de la red reducida.

get_keywords(reduce_graph(g, cutoff_cum))


# Puntos de cambio ----
cutoff_change <- find_cutoff(g, method="changepoint", knot_num=3)

cutoff_change

## Grafica
cutoff_fig +
  geom_hline(yintercept=cutoff_change, linetype="dashed")

## Eleccion de corte
g_redux <- reduce_graph(g, cutoff_change[1])
selected_terms <- get_keywords(g_redux)

selected_terms

## Terminos extras
extra_terms <- c(
  "medication",
  "cognitive-behavioural therapy",
  "cognitive behavioural therapy"
)

selected_terms <- c(selected_terms, extra_terms)

selected_terms


# Agrupamiento ----

grouped_terms <-list(
  medication=selected_terms[c(14, 28)],
  cbt=selected_terms[c(4, 6, 7, 17, 24, 25, 26, 29, 30)],
  phobia=selected_terms[c(2, 3, 9, 12, 15, 19, 20, 21, 23, 27)]
)

grouped_terms


write_search(
  grouped_terms,
  languages="English",
  exactphrase=TRUE,
  stemming=FALSE,
  closure="left",
  writesearch=TRUE
)

cat(read_file("search-inEnglish.txt"))

# Comprobando la nueva busqueda

new_results <- import_results(file="Datos/pubmed-panicORbeh-set.nbib")
nrow(new_results)


naive_results %>%
  dplyr::mutate(in_new_results=title %in% new_results[, "title"]) ->
  naive_results

naive_results %>%
  filter(!in_new_results) %>%
  select(title, keywords)

important_titles <- c(
  "Efficacy of treatments for anxiety disorder: A meta-analysis",
  "Cognitive behaviour therapy for health anxiety: A systematic review and meta-analysis",
  "A systematic review and meta-analysis of treatments for agrophobia"
)

data.frame(check_recall(important_titles, new_results[, "title"]))
