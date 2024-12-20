library(castor, quietly = TRUE)

Args <- commandArgs(TRUE)

full_tree <- read_tree(file=Args[1], check_label_uniqueness = TRUE)
known_tips <- read.table(Args[2], header=FALSE, stringsAsFactors = FALSE)$V1
output_path <- Args[3]

unknown_tips_index <- which(! full_tree$tip.label %in% known_tips)
unknown_tips <- full_tree$tip.label[unknown_tips_index]
all_tip_range <- 1:length(full_tree$tip.label)
known_tip_range <- which(! full_tree$tip.label %in% unknown_tips)

nsti_values <- find_nearest_tips(full_tree,
                                 target_tips=known_tip_range,
                                 check_input=TRUE)$nearest_distance_per_tip[unknown_tips_index]

nearest_tip <- find_nearest_tips(full_tree,
                                 target_tips=known_tip_range,
                                 check_input=TRUE)$nearest_tip_per_tip[unknown_tips_index]

nearest_tips_label <- full_tree$tip.label[nearest_tip]

write.table(x = data.frame("sequence" = unknown_tips, "NSTI" = nsti_values, "nearest_tip" = nearest_tips_label),
            file = output_path,
            sep="\t",
            quote = FALSE,
            row.names=FALSE)
