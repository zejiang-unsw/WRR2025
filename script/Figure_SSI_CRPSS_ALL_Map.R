# stationary model as the baseline model
rm(list=ls())
graphics.off() # clear environment and graphics

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

# Load required packages
#library(GPST)
library(WASP2.0)
library(hydroTSM) #time2season
library(sf)
library(ggplot2)
library(extRemes)
library(dplyr)
library(tidyr)
library(data.table)
library(export)
library(RColorBrewer)

library(scoringRules)

flag.save <- F

# Global parameters----
#display.brewer.all(colorblindFriendly = T, type=c("div","qual","seq","all")[2])
col_len <- 8
col_qual <- brewer.pal(n=col_len, name="Dark2")
col_pair <- brewer.pal(n=col_len+4, name="Paired")
col_div <- brewer.pal(n=col_len, name="RdYlBu")
col_seq <- brewer.pal(n=col_len, name="YlOrRd")

col_qual <- viridisLite::viridis(n=4, option = "D")

#labels_mod <- c("Stationary","Non-Stationary","Non-stationary (WASP)")
# labels_mod_exp <- c(expression(italic(mu~" ")), expression(italic(mu(x))),
#                     expression(italic(mu("x'"))), expression(italic(mu("x'(φ)"))))
labels_mod_exp <- c(expression(mu~" "), expression(mu(x)),
                    expression(mu("x'")), expression(mu("x'(φ)")))
lead_all <- c(1:3)
season_all <- c("DJF", "MAM", "JJA", "SON")

#flag.k <- 5; flag.iter= 3
flag.k <- 3; flag.iter= 1
log.tran <- 1
flag.v <- switch(1, paste0("_k",flag.k,"_log",log.tran), "_v1", "_v2","_v3") # k=5, iter=1

flag.v1 <- 1
flag.var <- switch(flag.v1, "wet","dry")
#flag.acc <- switch(flag.v1, 7, 30)

## WASP ----
mode <- switch(1, "DWT","MODWT") ## must use MODWT
if(mode=="DWT") wf <- switch(3, "haar","d4","d16") else wf <- "d4"

if(mode=="DWT") {method = "modwt"; pad = "sym"} # !!symmetric padding very important!!
boundary = switch(1,"periodic","reflection"); sign = "pos" # !! pos matters!!

## GEV -----
sig_lev <- 0.05 # likelihood ratio test
method_fevd <- switch(1, "MLE","GMLE") # !!! reliable parameter estimation is important !!!
#optim_args <- list(method=c("BFGS", "Nelder-Mead")[1])

# Load data---------------------------------------------------------------------
if(flag.var=="wet"){
  data("Qmax_SSI_AU_wet")
  Qmax_SSI_AU <- Qmax_SSI_AU_wet
} else {
  data("Qmax_SSI_AU_dry")
  Qmax_SSI_AU <- Qmax_SSI_AU_dry
}
data_q <- Qmax_SSI_AU 
stn_nos <- names(data_q)
CI_names <- c('Nino3', 'Nino4', 'Nino3.4', 'EMI', 'DMI','EPI', 'WPI', 'II', 
              'TSI', 'TBV', 'Disc')

CI_longname <-c('Niño 3', 'Niño 4', 'Niño 3.4', 'El Niño Modoki Index', 
                'Dipole Mode Index','Indian Ocean East Pole Index', 
                'Indian Ocean West Pole Index', 'Indonesian Index', 
                'Tasman Sea Index', 'Tropical Trans-basin Variability Index')

stn_check <- stn_nos
stn_samp <- sapply(stn_check, function(v) which(v==stn_nos))

## RPSS ----
file_save <- paste0("../result/SSI_CRPSS_",mode,"_",wf,"_Location_",
                    method_fevd,flag.v,"_ALL_AU_",flag.var,".Rdata")
load(file_save)

summary(RPSS_df_all)
RPSS_df_all1 <- RPSS_df_all %>% spread(Model, RPSS) %>% 
  mutate(Mod1=(Mod1-Mod0)/(1-Mod0)) %>% 
  mutate(Mod2=(Mod2-Mod0)/(1-Mod0)) %>% 
  mutate(Mod3=(Mod3-Mod0)/(1-Mod0)) %>% dplyr::select(!Mod0)

RPSS_df_all1$CI <- factor(RPSS_df_all1$CI, levels=CI_names)

# average across two folds
RPSS_df <- RPSS_df_all1 %>% dplyr::group_by(CI, Lead, Station) %>% 
  summarise(Mod1=mean(Mod1), Mod2=mean(Mod2), Mod3=mean(Mod3))
summary(RPSS_df)

RPSS_dt <- RPSS_df %>%  data.table()

##improvement -----
# lon and lat
filename <- "GRDC_Stations_Australia.csv"
excel_GRDC <- read.csv(paste0("../data/", filename),header = TRUE,
                       stringsAsFactors = FALSE)
ind_stn <- sapply(stn_nos, function(x) which(x== excel_GRDC[,1]))
#excel_GRDC[ind_stn,1] -  as.numeric(stn_nos)

loc_stn <- data.frame(Station=stn_nos, excel_GRDC[ind_stn,c("LONG_NEW","LAT_NEW")])

world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]

### map of predictor with averaged RPSS improvement for each station by each lead
summary(RPSS_dt)
RPSS_dt_pct1 <- NULL
p_map_list <- NULL
for(i_mod in 1:3){
  
  if(i_mod==1) {RPSS_dt1 <- RPSS_dt %>% subset(Mod2>Mod1 & Mod1>0) %>% mutate(Dif=(Mod2-Mod1))
  } else if(i_mod==2){
    RPSS_dt1 <- RPSS_dt %>% subset(Mod3>Mod1 & Mod1>0) %>% mutate(Dif=(Mod3-Mod1))
  } else {
    RPSS_dt1 <- RPSS_dt %>% subset(Mod3>Mod2 & Mod2>0) %>% mutate(Dif=(Mod3-Mod2))
  }

### maximum or mean improvement ----
RPSS_dt_avg <- RPSS_dt1[,.(RPSS=max(Dif)), by=c("Lead","Station")]
summary(RPSS_dt_avg)

RPSS_dt_sub <- RPSS_dt_avg %>% merge(loc_stn, by="Station") %>% data.table() %>% mutate(flag=1) %>%
  mutate(RPSS_cut = cut(RPSS, breaks= c(0, 0.05, 0.1, 0.2, Inf)))
labels_pct <- c("0~5%","5~10%","10~20%",">20%")

summary(RPSS_dt_sub)

## improvements by lead ----
tol_no <- 219
RPSS_dt_pct <- RPSS_dt_sub[,.(Sum=sum(flag)/tol_no),by=c("Lead")]
RPSS_dt_pct %>% arrange(Lead) %>% print()

#Map ----
loc_stn_sf <- sf::st_as_sf(x=RPSS_dt_sub, coords=c("LONG_NEW","LAT_NEW"), crs = st_crs(4326))
summary(loc_stn_sf)
loc_stn_sf$Lead <- factor(loc_stn_sf$Lead, labels = paste0("Lead ",0:2))

summary(RPSS_dt_sub)
p_map <- ggplot(data = loc_stn_sf) +
  #geom_sf(data = map, show.legend = FALSE) +
  #coord_sf(xlim=range(st_coordinates(map)[,1]), ylim=range(st_coordinates(map)[,2])) +
  
  #geom_point(aes(x=long, y=lat, color=RPSS_cut,size=RPSS_cut), shape=16, alpha=1) +
  
  geom_sf(aes(fill=RPSS_cut,size=RPSS_cut), shape=21, alpha=0.8, stroke = 0.1) +
  geom_sf(data=world, aes(group = continent), color = 'black', alpha=0.2) +
  #
  #geom_sf(data = lr_df_samp, aes(color=Station,shape=Station), size=3, stroke = 3, show.legend = TRUE) +
  facet_wrap(Lead~., ncol=1) +
  
  #coord_sf(xlim = c(-10,-60), ylim = c(-120,-170), expand = FALSE) +
  coord_sf(
    crs = st_crs(3577),  # Australian Albers projection
    xlim = c(-2050000, 2200000),  # Australia bounding box in EPSG 3577
    ylim = c(-5000000, -1000000),
    label_axes = list(
      bottom = "E", top = NA,
      left = "N", right = "N"
    )) +
  
  # scale_fill_manual(values = c(NA, NA, NA, NA,
  #                              NA, NA, NA, NA),
  #                   guide = 'none',
  #                   na.value = 'white') +
  #scale_size_continuous(limits = c(0,1), breaks = c(0,0.25,0.50, 0.75, 1)) +
  #scale_color_gradientn(colours=c("white",col_teal), limits=c(0,1), breaks = c(0,0.25,0.50, 0.75, 1)) +
  #scale_color_continuous(limits=c(0,1), breaks = c(0,0.25,0.50, 0.75, 1)) +
  scale_size_manual(values=c(0.5,1,1.5,2,2.5), label=labels_pct) +
  scale_fill_manual(values=col_pair[3:6], label=labels_pct) +
  
  labs(x=NULL, y=NULL, fill="Improvements in CRPSS", 
       size="Improvements in CRPSS") +
  
  egg::theme_article(base_size = 7) +
  theme(text = element_text(size = 7,family="sans"),
        plot.margin = unit(c(0.1,0.2,0.1,0.2), "cm"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linewidth = 0.05,
                                        color = alpha('grey',0.8)),
        strip.text.x = element_text(size = 7,family="sans", face='bold'),
        axis.ticks.length  = unit(0.1, "cm"),
        axis.ticks = element_line(linewidth = 0.2),
        
        #strip.background = element_blank(),
        strip.background = element_rect(fill="#DDDDDD", colour="grey4"),
        strip.placement = "outside",
        #strip.text.x = element_blank(),
        #strip.text.y = element_text(face='bold'),
        
        axis.text = element_text(size=7, family="sans", face="plain",color="black"),
        legend.position = "bottom",
        #legend.position = c(0.85,0.15),
        legend.direction = c("horizontal","vertical")[1],
        legend.key.spacing.x =  unit(0.1, "cm"),
        legend.text = element_text(size=8, margin = margin(l = -5)),
        legend.box.spacing = unit(0, "pt")
  )

p_map %>% print()

p_map_list[[i_mod]] <- p_map
### Figure a ----
if(flag.save){
filen4b <- paste0("Figure5b_CRPSS_max_m",i_mod,flag.v,"_AU_SSI_",flag.var,".png")
png(filen4b, height=7, width=18, units="cm", bg = "white", res=600)
p_map %>% print()
dev.off()

# graph2svg(x=p_map, file=filen, aspectr=2, font = "Arial", width = 9/2.54, # cm to inch
#           height = 7/2.54, bg = "white")
# graph2pdf(x=p_map, file=filen4b, aspectr=2, font = "Arial", width = 18/2.54,
#           height = 7/2.54, bg = "transparent")
# graph2eps(x=p_map, file=filen, aspectr=2, font = "Arial", width = 18/2.54,
#           height = 7/2.54, , bg = "transparent")
}

RPSS_dt_pct1 <- rbind(RPSS_dt_pct1, data.frame(Group=i_mod,
                                               RPSS_dt_sub[,.(pct=sum(flag)/tol_no),by=c("Lead","RPSS_cut")]))
#RPSS_dt_pct1 %>% arrange(Lead,RPSS_cut) %>% print()

}

### Figure b ----
## improvements by lead & group ----
summary(RPSS_dt_pct1)
RPSS_dt_pct1 %>% group_by(Lead,Group) %>% summarise(pct=sum(pct)) %>% spread(Group, pct) %>% print()

RPSS_dt_pct2 <- RPSS_dt_pct1 %>% mutate(RPSS=as.numeric(RPSS_cut))
RPSS_dt_pct2$RPSS_cut <- factor(RPSS_dt_pct2$RPSS_cut, labels=1:4)
summary(RPSS_dt_pct2)
head(RPSS_dt_pct2); head(RPSS_dt_pct1)

RPSS_dt_pct2 %>% subset(as.numeric(RPSS_cut)>2) %>% group_by(Lead,Group) %>% 
  summarise(pct=sum(pct)) %>% spread(Group, pct) %>% print()

RPSS_dt_pct1$Group <- factor(RPSS_dt_pct1$Group)
RPSS_dt_pct1$Lead <- factor(RPSS_dt_pct1$Lead, labels=paste0("Lead ",0:2))

p_bar <- ggplot(data=RPSS_dt_pct1 %>% subset(Group==2), 
                aes(x=RPSS_cut, y=pct*100, fill=RPSS_cut)) +
  
  geom_bar(position="stack", stat="identity", width = 0.7) + 

  #scale_color_manual(values=col_qual[c(3,2,1)], labels=labels_mod_exp[-1]) +
  scale_fill_manual(values=col_pair[3:6], labels=labels_pct) + 
  
  facet_wrap(.~Lead, scales ="fixed", ncol=1) +
  #facet_grid(RPSS_cut~Lead, scales ="free") +
  
  #scale_x_discrete(labels=labels_mod_exp[3:4]) +
  scale_x_discrete(labels=labels_pct) +
  scale_y_continuous(limits=c(0,30), breaks = scales::pretty_breaks(n = 5), 
                     expand=c(0,0)) +
  labs(x=NULL,y="Percentage of stations shown improvements (%)", 
       fill="Improvements \n in CRPSS (%)") +
  # guides(fill = guide_legend(title = "Lead",
  #                            label.position = "bottom",
  #                            title.position = "left", title.vjust = 1)) +
  
  guides(color = guide_legend(reverse = TRUE)) +
  
  egg::theme_article() +
  theme(text = element_text(size = 7,family="sans"),
        plot.margin = unit(c(0.1,0.5,0.1,0.5), "cm"),
        
        #panel.spacing.y = unit(1, "lines"), # space between subplot
        #strip.background = element_blank(),
        strip.background = element_rect(fill="#DDDDDD", colour="grey4"),
        strip.placement = "outside",
        strip.text.x = element_text(size = 7,family="sans", face='bold'),
        #strip.text.y = element_text(face='bold'),
        #strip.text.x.top = element_blank(),
        axis.text = element_text(size=7, family="sans", face="plain",color="black"),
        axis.title= element_text(size=7, family="sans", face="plain"),
        
        axis.title.x = element_text(margin = margin(t = 5)), 
        
        axis.ticks.length  = unit(0.1, "cm"),
        axis.ticks = element_line(linewidth = 0.2),
        
        legend.position = "none",
        #legend.position = c(-0.35,0.15),
        legend.title = element_text(hjust = 0.8),
        legend.direction = c("horizontal","vertical")[2],
        legend.background = element_rect(fill="transparent"),
        legend.key.spacing.x =  unit(0.1, "cm"),
        legend.key.spacing.y =  unit(0.1, "cm"),
        legend.key.width =  unit(0.5, "cm"),
        legend.text = element_text(size=7, margin = margin(r = 1, l=5))
        
  )
p_bar %>% print()

if(flag.save){
filen5c <- paste0("Figure5c_CRPSS_pct_m",flag.v,"_AU_SSI_",flag.var,".png")
png(filen5c, height=8, width=12, units="cm", bg = "white", res=600)
p_bar %>% print()
dev.off()
}


fig <- cowplot::plot_grid(p_map_list[[2]], p_bar, ncol=2,
                          labels = letters,
                          label_size = 10)
fig %>% print()

#shell.exec(filen5c)

# Lead     `1`   `2`   `3`
# <fct>  <dbl> <dbl> <dbl>
# 1 Lead 0 0.799 0.694 0.671
# 2 Lead 1 0.767 0.626 0.525
# 3 Lead 2 0.822 0.566 0.329
