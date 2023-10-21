library(dplyr)
library(stringr)

file_name <- 'result_guide2'
offtargets <- read.csv(paste (file_name, '.txt', sep = ''), sep = '\t')

#ищем оффтаргеты с только одним мисматчем. Если всего один мисматч, то проверить, не находится ли он в ближайших четырех нуклеотидах от ПАМ. Если он не там, то ставим 1
offtargets$onemismatch <- ifelse (grepl("^[A-Z]*[a-z\\-]{1}[A-Z]*$", offtargets$DNA) , 
                                  ifelse (grepl("^[^a-z]*[a-z\\-]{1}[^a-z]*$", str_sub(offtargets$DNA, -7)), 0, 1), 0)
#ищем оффтаргеты, где нет мисматчей или пропусков в 8членном коре. Если такие есть, проверяем, есть ли у них DNA bulge в 8 членном коре. Если нет, ставим 1
offtargets$no8mismatch <- ifelse (!grepl("[a-z\\-]",str_sub(offtargets$DNA, -11)),
                                   ifelse (!grepl("[-]",str_sub(offtargets$crRNA, -11)), 1, 0), 0)

offtargets <- offtargets[order(-offtargets$no8mismatch, -offtargets$onemismatch, -offtargets$Position),]
#ищем оффтаргеты, где сумма мисматчей и балджей - меньше или равна 5.
offtargets$sum <- ifelse (offtargets$Mismatches + offtargets$Bulge.Size <=5, 1, 0)
#Важные оффтаргеты - это те, у которых сумма мисматчей и пропусков меньше или равна 5 И (либо нет мисматчей в 8 членном коре, либо один мисматч за пределами 4 членного кора)
offtargets$isimportant <- ifelse (offtargets$sum==1&(offtargets$no8mismatch == 1|offtargets$onemismatch == 1), 1, 0)
offtargets$Short_Position <- as.numeric (str_extract(paste(offtargets$Chromosome, offtargets$Position, sep = ''), '[0-9]+(?=.)'))
offtargets <- offtargets[order(-offtargets$isimportant, -offtargets$Short_Position),]
important_offtargets <- subset (offtargets[!duplicated(offtargets[,'Short_Position']),], isimportant == 1)
important_offtargets$for_blast <- str_remove_all (paste ('>', c(1:nrow(important_offtargets)), '\n', important_offtargets$DNA, sep = ''), "-") 

short  <- select (important_offtargets, for_blast)

write.csv(short, paste ('dangerous', file_name, '.txt'), row.names = FALSE, quote=FALSE)
