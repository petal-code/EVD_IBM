library(dplyr)
library(ggplot2)

inf_df <- result$infected
nodes  <- readRDS(sprintf("output/network/%s_nodes.rds", case_tag))

# Household size per person
hh_size_df <- nodes %>%
  group_by(hh_id) %>%
  mutate(hh_size = n()) %>%
  ungroup() %>%
  select(person_id, hh_id, hh_size)

# For each infected person, find HH contacts and SAR
layer1 <- readRDS(sprintf("output/network/%s_layer1_household.rds", case_tag))

hh_members <- bind_rows(
  layer1 %>% select(ego=from, alter=to),
  layer1 %>% select(ego=to,   alter=from)
)

# Infected set
infected_ids <- inf_df$person_id

# Per household: index case (HH route) → SAR among HH members
hh_infected <- inf_df %>%
  filter(contact_type == 1L) %>%   # HH transmission
  left_join(hh_size_df, by="person_id")

# For each HH-infected person, what fraction of their HH got infected?
sar_df <- hh_members %>%
  filter(ego %in% infected_ids) %>%
  left_join(hh_size_df %>% select(person_id, hh_size),
            by=c("ego"="person_id")) %>%
  group_by(ego, hh_size) %>%
  summarise(
    n_hh_contacts  = n(),
    n_hh_infected  = sum(alter %in% infected_ids),
    SAR            = n_hh_infected / n_hh_contacts,
    .groups        = "drop"
  ) %>%
  filter(n_hh_contacts > 0)

# Plot
p_sar <- ggplot(sar_df, aes(x=factor(hh_size), y=SAR)) +
  geom_jitter(width=0.2, height=0.01, alpha=0.3, size=0.8, color="steelblue") +
  stat_summary(fun=mean, geom="point", color="red", size=3) +
  stat_summary(fun=mean, geom="line", aes(group=1), color="red", linewidth=0.8) +
  scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
  labs(title    = sprintf("%s — SAR by household size (w=%.1f)", case_tag, w_household),
       subtitle = sprintf("Red = mean SAR | N infected = %d", nrow(inf_df)),
       x="Household size", y="Secondary Attack Rate") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"))

print(p_sar)
ggsave(sprintf("figure/C1_sim/%s_SAR_by_hhsize.png", case_tag),
       p_sar, width=8, height=5, dpi=150)
