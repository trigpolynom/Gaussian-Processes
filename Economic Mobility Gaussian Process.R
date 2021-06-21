library(ggplot2)
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)
library(stringi)
library(tidyverse)
library(GGally)
library(earth)
library(vip)
library(pdp)
library(pls)
library(visdat)
library(gbm)
library(tmap)
library(sf)
library(spdep)
library(rgdal)
library(rmapshaper)
library(tidycensus)
library(gstat)
library(INLA)
library(viridis)
library(grid)
library(maps)
library(ggthemes)
library(scales)
library(dismo)
library(fitdistrplus)
library(maptools)
library(ggmap)
library(matrixStats)
library(geoR)
library(automap)
library(stars)
library(visdat)
library(dplyr)
library(raster)
library(lattice)
library(spatialreg)


dss_index <- st_as_sf(Dataset_for_Ethan,
                      coords = c("Longitude_POINT_X", "Lattitude_POINT_Y"),
                      crs = 27561)

summary(dss_index)

plot(dss_index$geometry)

rn <- row.names(Dataset_for_Ethan)
k1 <- knn2nb(knearneigh(dss_coords))
all.linked <- max(unlist(nbdists(k1, dss_coords)))
col.nb.0.all <- dnearneigh(dss_coords, 0, all.linked, row.names=rn)
DSS.W <- nb2listw(col.nb.0.all, style = "W", zero.policy = T) 

summary(col.nb.0.all, dss_coords)

summary(DSS.W, zero.policy=TRUE)


dss_OLS <- lm(d_Diss_Index~d_Municipalities+
                d_share_black+
                d_Share_Urban+
                d_Pop_Growth+
                d_Poor+
                d_Thiel_Index_MedHHInc+
                d_Vacant_Houses+
                d_TotCrimeRate+
                d_HomeValue+
                d_FundsPerStudent+
                d_DrugMortRate+
                d_EmpRate+
                d_Renters+
                d_PublicTrans,
              data = dss_index)

summary(dss_OLS)

##Moran Test

moran <- lm.morantest(model = dss_OLS, 
                      listw = DSS.W,
                      zero.policy = TRUE)

print(moran)

##Lagrange Multiplier Diagnostics

LMtest1 <- lm.LMtests(dss_OLS,
                      DSS.W,
                      zero.policy = TRUE,
                      test = c('LMlag', 'LMerr'))

print(LMtest1)

LMtest2 <- lm.LMtests(dss_OLS,
                      DSS.W,
                      zero.policy = TRUE,
                      test = c('RLMlag', 'RLMerr'))
print(LMtest2)

##Spatial Error Model
dss_err <- errorsarlm(d_Diss_Index ~ d_Municipalities+
                        d_share_black+
                        d_Share_Urban+
                        d_Pop_Growth+
                        d_Poor+
                        d_Thiel_Index_MedHHInc+
                        d_Vacant_Houses+
                        d_TotCrimeRate+
                        d_HomeValue+
                        d_FundsPerStudent+
                        d_DrugMortRate+
                        d_EmpRate+
                        d_Renters+
                        d_PublicTrans,
                      data = dss_index,
                      DSS.W,
                      zero.policy = TRUE)

summary(dss_err)

##Spatial Lag Model

dss_lag <- lagsarlm(d_Diss_Index ~ d_Municipalities+
                      d_share_black+
                      d_Share_Urban+
                      d_Pop_Growth+
                      d_Poor+
                      d_Thiel_Index_MedHHInc+
                      d_Vacant_Houses+
                      d_TotCrimeRate+
                      d_HomeValue+
                      d_FundsPerStudent+
                      d_DrugMortRate+
                      d_EmpRate+
                      d_Renters+
                      d_PublicTrans,
                    data = dss_index,
                    DSS.W,
                    zero.policy = TRUE)

summary(dss_lag)

##impacts

lag.impacts <- impacts(dss_lag, listw = DSS.W)
print(lag.impacts)

plot(st_geometry(dss_index), border="grey", reset=FALSE,
     main=paste("Distance based neighbours 0-",  format(all.linked), sep=""))

plot(col.nb.0.all, dss_coords, add=TRUE)

tmap_mode('view')
tmap_options(check.and.fix = TRUE)
tm_shape(dss_index) + tm_dots(col = "d_Diss_Index", size = 0.05, midpoint = NA)

#spatial weights matrix
dss_coords <- st_coordinates(dss_index)
dss.nb <- dnearneigh(dss_coords,
                         0,
                     50000,
                         longlat = TRUE)


summary(dss.W, zero.policy=TRUE)

if(requireNamespace("gstat", quietly = TRUE)){
  data(coalash, package="gstat")
   summary(coalash)
  } else { tmap_options(check.and.fix = TRUE) and rerun the plot
     install.packages("gstat")
     requireNamespace("gstat", quietly = TRUE)
     data(coalash, package="gstat")
     summary(coalash)
  }

##NC##

coordinates(For_Ethan_LattsANDLongs) <- c("Longitude_POINT_X", "Lattitude_POINT_Y")
coordinates(coords_test) <- c("POINT_X_longitude", "POINT_Y_latitude")

nc_full <- For_Ethan_LattsANDLongs %>%
  drop_na()

e <- extent(bbox(For_Ethan_LattsANDLongs))

e <- e + 1000

r <- raster(e)

x <- rasterize()

plot(r)

initialsplitkrige <- initial_split(as.data.frame(For_Ethan_LattsANDLongs))
nc_train <- training(initialsplitkrige)
nc_test <- testing(initialsplitkrige)

var1 <- autofitVariogram(kir_top20_pooled_pooled_p25~1, nc_train)

kriging1 <- autoKrige(kir_top20_pooled_pooled_p25~1,
                      nc_train)

kriging2 <- autoKrige(kir_top20_pooled_pooled_p25~1,
                      nc_train,
                      nc_test,
                      fix.values = c(0.0015,NA,NA))

kriging1 <- autoKrige.cv(kir_top20_pooled_pooled_p25~1,
                      nc_train)

kriging2 <- autoKrige.cv(kir_top20_pooled_pooled_p25~Share_frac_coll_plus2000+,
                         nc_train,
                         nfold = 10)

summary.autoKrige.cv(kriging2)

automapPlot(kriging1$krige_output, "var1.pred", 
            sp.layout = list("sp.points", nc_train))

compare.cv()
summary.autoKrige.cv = function(object, ..., digits = 4) {
  # The summary function for the autoKrige.cv object. Returns some
  # statistics about the cross-validation
  obj = object$krige.cv_output
  out = list()
  # mean error, ideally 0:
  out$mean_error = mean(obj$residual)
  # mean error divided by the mean of the observed values, measure for how large the mean_error is in contrast to the mean of the dataset
  out$me_mean = out$mean_error / mean(obj$observed)
  # mean absolute error, ideally 0, less vulnerable to outliers
  out$MAE = mean(abs(obj$residual))
  # MSE, ideally small
  out$MSE = mean(obj$residual^2)
  # Mean square normalized error, ideally close to 1
  out$MSNE = mean(obj$zscore^2)
  # correlation observed and predicted, ideally 1
  out$cor_obspred = cor(obj$observed, obj$observed - obj$residual)
  # correlation predicted and residual, ideally 0
  out$cor_predres = cor(obj$observed - obj$residual, obj$residual)
  # RMSE, ideally small
  out$RMSE = sqrt(sum(obj$residual^2) / length(obj$residual))
  # RMSE / sd(observed), measure for how much the residuals vary to the total variation in the dataset
  out$RMSE_sd = out$RMSE / sd(obj$observed)
  # URMSE, ideally zero
  out$URMSE = sqrt((sum(obj$residual^2) / length(obj$residual)) - mean(obj$residual)^2)
  # Inter quartile range, ideally small
  out$iqr = IQR(obj$residual)
  
  out = lapply(out, signif, digits = digits)
  out = t(t(out))
  return(out)
}


models1 <- lm(kir_top20_pooled_pooled_p25 ~ hhinc_mean2000+
                mean_commutetime2000+
                Share_frac_coll_plus2000+
                Share_foreign_share2010+
                popdensity2000+
                Share_Poor_2000+
                Share_Black_2000+
                Share_SingleParent_2000+
                job_density_2013+
                Share_Commutes_Less15mins_2010+
                Share_Commuters_Over1Hr_2010+
                Share_Commuters_PubTrans_2010+
                Share_Commuters_PersVeh_2010+
                Share_HH_withPub_Assist_Inc_2010+
                Share_Industry_AgForFishMin_2010+
                Share_Industry_Construction_2010+
                Share_Industry_Manufact_2010+
                Share_Industry_RetailTrade_2010+
                Share_Industry_Information_2010+
                Share_Industry_Finance_2010+
                Share_Indus_EdcSrvHelSocAs_2010+
                Share_Industry_WhlesleTrade_2010+
                Share_Industry_TrnsWareUtil_2010+
                Share_Industry_PrfSciManWst_2010+
                Share_Indus_ArtEntrRecAcc_2010+
                Share_Indus_PubAdmin_2010+
                Share_Vacant_Units_2010+
                Gini_Index_2010+
                Male_life_expect_2010+
                Share_Heavy_Drinker_2010+
                Total_Lynchings_1983to1941+
                Entropy_Score_2000+
                Share_Asian_2000_IMP+
                rent_burden_2010_IMP+
                Tract_Had_Redlines +
                Med_Home_Value_2010_IMP+
                Ln_Dist_NCAandT_Univ+
                Ln_NC_CENTRAL_UNIV+
                Ln_Dist_WSSU+
                Share_Dem_Vote_2000+
                Unemp_Rate_2000+
                Real_GDP_2012+
                Share_Homeownership_2000+
                gsmn_math_g3_2013+
                UNIV_COUNT+
                Net_Mig_Rate_NWs_1970s+
                Net_Mig_Rate_NWs_1980s, 
              data = For_Ethan_LattsANDLongs@data)
summary(models1)
hist(residuals(models1))
rug(residuals(models1))
ix <- which.max(residuals(models1))
For_Ethan_LattsANDLongs@data[ix,]

RSS <- c(crossprod(models1$residuals))

MSE <- RSS / length(models1$residuals)

RMSE <- sqrt(MSE)

summary.lm(models1)

neighbor_list <- nb2listw(neib, style = "B", zero.policy = TRUE)

(morans <- lm.morantest(models1, neighbor_list, zero.policy = TRUE))

w2 <- knn2nb(knearneigh(For_Ethan_LattsANDLongs, k=8))
moran.test(For_Ethan_LattsANDLongs$kir_top20_pooled_pooled_p25, nb2listw(w2))

n <- 7
res <- data.frame(k=2^(1:n), I=rep(NA,n))

for(i in 1:n){
  w <- knn2nb(knearneigh(For_Ethan_LattsANDLongs, k=2^i))
  res$I[i] <- moran.test(For_Ethan_LattsANDLongs$kir_top20_pooled_pooled_p25, nb2listw(w))$estimate[1]
}

plot(res, type="b",
     main="Moran's I by number of neighbors",
     pch= 20,
     cex= 1.5)

library(ncf)
mobilityI <- spline.correlog(x= coordinates(For_Ethan_LattsANDLongs)[,],
                         y= coordinates(For_Ethan_LattsANDLongs)[,2],
                         z= For_Ethan_LattsANDLongs$kir_top20_pooled_pooled_p25,
                         resamp= 50,
                         xmax= 10000,
                         quiet= TRUE)

plot(mobilityI)

dim(For_Ethan_LattsANDLongs)

ac_nb2 <- poly2nb(AC, queen = FALSE)

class(For_Ethan_LattsANDLongs)

proj4string(For_Ethan_LattsANDLongs) <- CRS("+init=epsg:4326")

is.projected(For_Ethan_LattsANDLongs)

plot(For_Ethan_LattsANDLongs, border="grey60", axes=TRUE, asp=1)

mobility_sf <- st_as_sf(For_Ethan_LattsANDLongs,
                        coords = c("Longitude_POINT_X", "Lattitude_POINT_Y"))

writeOGR(For_Ethan_LattsANDLongs, "CS/ML/PhillyHomicides", "PhillyHomcides", driver = "ESRI Shapefile", overwrite_layer = TRUE)

st_crs(mobility_sf) <- 4326

st_crs(mobility_sf)

class(mobility_sf)

mobility_sf.sf <- st_read(mobility_sf)

st_write(mobility_sf, "C:/Users/egordon/Desktop/CS/ML", driver = "ESRI Shapefile")

mobility <- as.data.frame(For_Ethan_LattsANDLongs)
Fina_Dataset_toUse %>%
  ggplot()+
  #geom_histogram(aes(kir_top20_pooled_pooled_p25), binwidth = 0.005, fill = "black", alpha = 0.2)+
  geom_histogram(aes(kir_top20_black_pooled_p25), binwidth = 0.005, fill = "red", alpha = 0.2)+
  geom_histogram(aes(kir_top20_hisp_pooled_p25), binwidth = 0.005, fill = "blue", alpha = 0.2)+
  geom_histogram(aes(kir_top20_white_pooled_p25), binwidth = 0.005, fill = "black", alpha = 0.3)

#Poverty
Fina_Dataset_toUse %>%
  select(kir_top20_pooled_pooled_p25,
         Share_Poor_2010,
         Share_Poor_2000,
         Share_Poor_1990) %>%
  ggpairs()

#Race
Fina_Dataset_toUse %>%
  select(kir_top20_pooled_pooled_p25,
         Share_Black_2010,
         Share_Black_2000,
         Share_white_2010,
         Share_White_2000,
         Share_Asian_2000_IMP,
         Share_Asian_2010_IMP) %>%
  ggpairs()

#SEC
Fina_Dataset_toUse %>%
  select(kir_top20_pooled_pooled_p25,
         Share_Black_2000,
         Share_frac_coll_plus2000,
         Share_frac_coll_plus2010,
         Share_SingleParent_1990,
         Share_SingleParent_2000,
         Share_SingleParent_2010) %>%
  pairs(lower.panel = panel.smooth, col = "blue", cex = 0.8)
  
##missings
Fina_Dataset_toUse %>%
  is.na() %>%
  reshape2::melt() %>%
  ggplot(aes(Var2, Var1, fill=value)) + 
  geom_raster() + 
  coord_flip() +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_grey(name = "", 
                  labels = c("Present", 
                             "Missing")) +
  xlab("Observation") +
  theme(axis.text.y  = element_text(size = 4))

vis_miss(Fina_Dataset_toUse, cluster = TRUE)

#nzv
nzv <- nearZeroVar(em_train, saveMetrics= TRUE)

nzv[nzv$nzv,][1:10,]


#cor
em_cor <- em_train %>%
  drop_na() %>%
  select_if(is.numeric) %>%
  cor()

view(Fina_Dataset_toUse %>%
  drop_na()%>%
  select_if(is.numeric) %>%
  cor())

summary(em_cor[upper.tri(em_cor)])

highly_corr <- em_train %>%
  drop_na() %>%
  select_if(is.numeric) %>%
  findCorrelation(cutoff = 0.75)

em_train1 <- em_train[,-highly_corr]

cor2 <- cor(em_train1)

summary(cor2[upper.tri(cor2)])

#comboInfo
comboinfo <- em_train %>%
  drop_na() %>%
  select_if(is.numeric) %>%
  findLinearCombos()

em_combo <- em_train[, -comboinfo$remove]

#principal component analysis

apply(Fina_Dataset_toUse, 2, var)

scaled_df <- Fina_Dataset_toUse %>%
  select_if(is.numeric) %>%
  apply(2, scale)

pca_result <- Fina_Dataset_toUse %>%
  drop_na() %>%
  select_if(is.numeric) %>%
  prcomp(scale = FALSE)

names(pca_result)
pca_result$center
pca_result$rotation
pca_result$x

biplot(pca_result, scale = 1)

#real pre process
preProcValues <- preProcess(em_train, method = c("center", "scale"))

trainTransformed <- predict(preProcValues, em_train)
testTransformed <- predict(preProcValues, em_test)


set.seed(123)  # for reproducibility
(cv_model1 <- train(
  form = kir_top20_pooled_pooled_p25 ~ hhinc_mean2000+
    mean_commutetime2000+
    Share_frac_coll_plus2000+
    Share_foreign_share2010+
    popdensity2000+
    Share_Poor_2000+
    Share_Black_2000+
    Share_SingleParent_2000+
    job_density_2013+
    Share_Commutes_Less15mins_2010+
    Share_Commuters_Over1Hr_2010+
    Share_Commuters_PubTrans_2010+
    Share_Commuters_PersVeh_2010+
    Share_HH_withPub_Assist_Inc_2010+
    Share_Industry_AgForFishMin_2010+
    Share_Industry_Construction_2010+
    Share_Industry_Manufact_2010+
    Share_Industry_RetailTrade_2010+
    Share_Industry_Information_2010+
    Share_Industry_Finance_2010+
    Share_Indus_EdcSrvHelSocAs_2010+
    Share_Industry_WhlesleTrade_2010+
    Share_Industry_TrnsWareUtil_2010+
    Share_Industry_PrfSciManWst_2010+
    Share_Indus_ArtEntrRecAcc_2010+
    Share_Indus_PubAdmin_2010+
    Share_Vacant_Units_2010+
    Gini_Index_2010+
    Male_life_expect_2010+
    Share_Heavy_Drinker_2010+
    Total_Lynchings_1983to1941+
    Entropy_Score_2000+
    `Fe,ale_obesity_prevalence_2009`+
    Share_Asian_2000_IMP+
    rent_burden_2010_IMP+
    Tract_Had_Redlines +
    Med_Home_Value_2010_IMP+
    Ln_Dist_NCAandT_Univ+
    Ln_NC_CENTRAL_UNIV+
    Ln_Dist_WSSU+
    Share_Dem_Vote_2000+
    Unemp_Rate_2000+
    Real_GDP_2012+
    Share_Homeownership_2000+
    gsmn_math_g3_2013+
    UNIV_COUNT+
    Net_Mig_Rate_NWs_1970s+
    Net_Mig_Rate_NWs_1980s, 
  data = em_train, 
  method = "lm",
  trControl = trainControl(method = "cv", number = 10)
))

#MARS
hyper_grid <- expand.grid(
  degree = 1:3, 
  nprune = seq(2, 100, length.out = 10) %>% floor()
)

hyper_grid1 <- expand.grid(
  degree = 1:4, 
  nprune = seq(6, 50, length.out = 10) %>% floor()
)

set.seed(123)

(tuned_mars <- train(kir_top20_pooled_pooled_p25 ~ hhinc_mean2000+
                    mean_commutetime2000+
                    Share_frac_coll_plus2000+
                    Share_foreign_share2010+
                    popdensity2000+
                    Share_Poor_2000+
                    Share_Black_2000+
                    Share_SingleParent_2000+
                    job_density_2013+
                    Share_Commutes_Less15mins_2010+
                    Share_Commuters_Over1Hr_2010+
                    Share_Commuters_PubTrans_2010+
                    Share_Commuters_PersVeh_2010+
                    Share_HH_withPub_Assist_Inc_2010+
                    Share_Industry_AgForFishMin_2010+
                    Share_Industry_Construction_2010+
                    Share_Industry_Manufact_2010+
                    Share_Industry_RetailTrade_2010+
                    Share_Industry_Information_2010+
                    Share_Industry_Finance_2010+
                    Share_Indus_EdcSrvHelSocAs_2010+
                    Share_Industry_WhlesleTrade_2010+
                    Share_Industry_TrnsWareUtil_2010+
                    Share_Industry_PrfSciManWst_2010+
                    Share_Indus_ArtEntrRecAcc_2010+
                    Share_Indus_PubAdmin_2010+
                    Share_Vacant_Units_2010+
                    Gini_Index_2010+
                    Male_life_expect_2010+
                    Share_Heavy_Drinker_2010+
                    Total_Lynchings_1983to1941+
                    Entropy_Score_2000+
                    `Fe,ale_obesity_prevalence_2009`+
                    Share_Asian_2000_IMP+
                    rent_burden_2010_IMP+
                    Tract_Had_Redlines +
                    Med_Home_Value_2010_IMP+
                    Ln_Dist_NCAandT_Univ+
                    Ln_NC_CENTRAL_UNIV+
                    Ln_Dist_WSSU+
                    Share_Dem_Vote_2000+
                    Unemp_Rate_2000+
                    Real_GDP_2012+
                    Share_Homeownership_2000+
                    gsmn_math_g3_2013+
                    UNIV_COUNT+
                    Net_Mig_Rate_NWs_1970s+
                    Net_Mig_Rate_NWs_1980s,
  data = em_train,
  method = "earth",
  metric = "RMSE",
  trControl = trainControl(method = "cv", number = 10),
  tuneGrid = hyper_grid
))


#black mobility

black_mobility <- Fina_Dataset_toUse %>%
  select(-c(kir_top20_pooled_pooled_p25,
            kir_top20_white_pooled_p25,
            kir_top20_hisp_pooled_p25)) %>%
  drop_na()


set.seed(123)
(tuned_black <- train(kir_top20_black_pooled_p25 ~ hhinc_mean2000+
                       mean_commutetime2000+
                       Share_frac_coll_plus2000+
                       Share_foreign_share2010+
                       popdensity2000+
                       Share_Poor_2000+
                       Share_Black_2000+
                       Share_SingleParent_2000+
                       job_density_2013+
                       Share_Commutes_Less15mins_2010+
                       Share_Commuters_Over1Hr_2010+
                       Share_Commuters_PubTrans_2010+
                       Share_Commuters_PersVeh_2010+
                       Share_HH_withPub_Assist_Inc_2010+
                       Share_Industry_AgForFishMin_2010+
                       Share_Industry_Construction_2010+
                       Share_Industry_Manufact_2010+
                       Share_Industry_RetailTrade_2010+
                       Share_Industry_Information_2010+
                       Share_Industry_Finance_2010+
                       Share_Indus_EdcSrvHelSocAs_2010+
                       Share_Industry_WhlesleTrade_2010+
                       Share_Industry_TrnsWareUtil_2010+
                       Share_Industry_PrfSciManWst_2010+
                       Share_Indus_ArtEntrRecAcc_2010+
                       Share_Indus_PubAdmin_2010+
                       Share_Vacant_Units_2010+
                       Gini_Index_2010+
                       Male_life_expect_2010+
                       Share_Heavy_Drinker_2010+
                       Total_Lynchings_1983to1941+
                       Entropy_Score_2000+
                       `Fe,ale_obesity_prevalence_2009`+
                       Share_Asian_2000_IMP+
                       rent_burden_2010_IMP+
                       Tract_Had_Redlines +
                       Med_Home_Value_2010_IMP+
                       Ln_Dist_NCAandT_Univ+
                       Ln_NC_CENTRAL_UNIV+
                       Ln_Dist_WSSU+
                       Share_Dem_Vote_2000+
                       Unemp_Rate_2000+
                       Real_GDP_2012+
                       Share_Homeownership_2000+
                       gsmn_math_g3_2013+
                       UNIV_COUNT+
                       Net_Mig_Rate_NWs_1970s+
                       Net_Mig_Rate_NWs_1980s,
                     data = black_mobility,
                     method = "earth",
                     metric = "RMSE",
                     trControl = trainControl(method = "cv", number = 10),
                     tuneGrid = hyper_grid
))


tuned_black$results %>%
  filter(nprune == tuned_black$bestTune$nprune, degree == tuned_black$bestTune$degree)


set.seed(123)
(cv_model1black <- train(kir_top20_black_pooled_p25 ~ hhinc_mean2000+
                     mean_commutetime2000+
                     Share_frac_coll_plus2000+
                     Share_foreign_share2010+
                     popdensity2000+
                     Share_Poor_2000+
                     Share_Black_2000+
                     Share_SingleParent_2000+
                     job_density_2013+
                     Share_Commutes_Less15mins_2010+
                     Share_Commuters_Over1Hr_2010+
                     Share_Commuters_PubTrans_2010+
                     Share_Commuters_PersVeh_2010+
                     Share_HH_withPub_Assist_Inc_2010+
                     Share_Industry_AgForFishMin_2010+
                     Share_Industry_Construction_2010+
                     Share_Industry_Manufact_2010+
                     Share_Industry_RetailTrade_2010+
                     Share_Industry_Information_2010+
                     Share_Industry_Finance_2010+
                     Share_Indus_EdcSrvHelSocAs_2010+
                     Share_Industry_WhlesleTrade_2010+
                     Share_Industry_TrnsWareUtil_2010+
                     Share_Industry_PrfSciManWst_2010+
                     Share_Indus_ArtEntrRecAcc_2010+
                     Share_Indus_PubAdmin_2010+
                     Share_Vacant_Units_2010+
                     Gini_Index_2010+
                     Male_life_expect_2010+
                     Share_Heavy_Drinker_2010+
                     Total_Lynchings_1983to1941+
                     Entropy_Score_2000+
                     `Fe,ale_obesity_prevalence_2009`+
                     Share_Asian_2000_IMP+
                     rent_burden_2010_IMP+
                     Tract_Had_Redlines +
                     Med_Home_Value_2010_IMP+
                     Ln_Dist_NCAandT_Univ+
                     Ln_NC_CENTRAL_UNIV+
                     Ln_Dist_WSSU+
                     Share_Dem_Vote_2000+
                     Unemp_Rate_2000+
                     Real_GDP_2012+
                     Share_Homeownership_2000+
                     gsmn_math_g3_2013+
                     UNIV_COUNT+
                     Net_Mig_Rate_NWs_1970s+
                     Net_Mig_Rate_NWs_1980s, 
                   data = black_mobility, 
                   method = "lm",
                   metric = "RMSE",
                   trControl = trainControl(method = "cv", number = 10)#,
                   #preProcess = c("zv", "center", "scale")
))


tuned_mars$results %>%
  filter(nprune == tuned_mars$bestTune$nprune, degree == tuned_mars$bestTune$degree)

cv_model1$results %>%
  filter(nprune == cv_model1$bestTune$nprune, degree == cv_model1$bestTune$degree)

#(black mobility) extract out of sample performance measures
summary(resamples(list(
  Multiple_regression = cv_model1black, 
  MARS = tuned_bm
)))$statistics$RMSE %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))


##white mobility

white_mobility <- Fina_Dataset_toUse %>%
  select(-c(kir_top20_pooled_pooled_p25,
            kir_top20_black_pooled_p25,
            kir_top20_hisp_pooled_p25)) %>%
  drop_na()


set.seed(123)
(tuned_white <- train(kir_top20_white_pooled_p25 ~ hhinc_mean2000+
                        mean_commutetime2000+
                        Share_frac_coll_plus2000+
                        Share_foreign_share2010+
                        popdensity2000+
                        Share_Poor_2000+
                        Share_Black_2000+
                        Share_SingleParent_2000+
                        job_density_2013+
                        Share_Commutes_Less15mins_2010+
                        Share_Commuters_Over1Hr_2010+
                        Share_Commuters_PubTrans_2010+
                        Share_Commuters_PersVeh_2010+
                        Share_HH_withPub_Assist_Inc_2010+
                        Share_Industry_AgForFishMin_2010+
                        Share_Industry_Construction_2010+
                        Share_Industry_Manufact_2010+
                        Share_Industry_RetailTrade_2010+
                        Share_Industry_Information_2010+
                        Share_Industry_Finance_2010+
                        Share_Indus_EdcSrvHelSocAs_2010+
                        Share_Industry_WhlesleTrade_2010+
                        Share_Industry_TrnsWareUtil_2010+
                        Share_Industry_PrfSciManWst_2010+
                        Share_Indus_ArtEntrRecAcc_2010+
                        Share_Indus_PubAdmin_2010+
                        Share_Vacant_Units_2010+
                        Gini_Index_2010+
                        Male_life_expect_2010+
                        Share_Heavy_Drinker_2010+
                        Total_Lynchings_1983to1941+
                        Entropy_Score_2000+
                        `Fe,ale_obesity_prevalence_2009`+
                        Share_Asian_2000_IMP+
                        rent_burden_2010_IMP+
                        Tract_Had_Redlines +
                        Med_Home_Value_2010_IMP+
                        Ln_Dist_NCAandT_Univ+
                        Ln_NC_CENTRAL_UNIV+
                        Ln_Dist_WSSU+
                        Share_Dem_Vote_2000+
                        Unemp_Rate_2000+
                        Real_GDP_2012+
                        Share_Homeownership_2000+
                        gsmn_math_g3_2013+
                        UNIV_COUNT+
                        Net_Mig_Rate_NWs_1970s+
                        Net_Mig_Rate_NWs_1980s,
                      data = white_mobility,
                      method = "earth",
                      metric = "RMSE",
                      trControl = trainControl(method = "cv", number = 10),
                      tuneGrid = hyper_grid
))


tuned_white$results %>%
  filter(nprune == tuned_bm$bestTune$nprune, degree == tuned_bm$bestTune$degree)


set.seed(123)
(cv_model1white <- train(kir_top20_white_pooled_p25 ~ hhinc_mean2000+
                           mean_commutetime2000+
                           Share_frac_coll_plus2000+
                           Share_foreign_share2010+
                           popdensity2000+
                           Share_Poor_2000+
                           Share_Black_2000+
                           Share_SingleParent_2000+
                           job_density_2013+
                           Share_Commutes_Less15mins_2010+
                           Share_Commuters_Over1Hr_2010+
                           Share_Commuters_PubTrans_2010+
                           Share_Commuters_PersVeh_2010+
                           Share_HH_withPub_Assist_Inc_2010+
                           Share_Industry_AgForFishMin_2010+
                           Share_Industry_Construction_2010+
                           Share_Industry_Manufact_2010+
                           Share_Industry_RetailTrade_2010+
                           Share_Industry_Information_2010+
                           Share_Industry_Finance_2010+
                           Share_Indus_EdcSrvHelSocAs_2010+
                           Share_Industry_WhlesleTrade_2010+
                           Share_Industry_TrnsWareUtil_2010+
                           Share_Industry_PrfSciManWst_2010+
                           Share_Indus_ArtEntrRecAcc_2010+
                           Share_Indus_PubAdmin_2010+
                           Share_Vacant_Units_2010+
                           Gini_Index_2010+
                           Male_life_expect_2010+
                           Share_Heavy_Drinker_2010+
                           Total_Lynchings_1983to1941+
                           Entropy_Score_2000+
                           `Fe,ale_obesity_prevalence_2009`+
                           Share_Asian_2000_IMP+
                           rent_burden_2010_IMP+
                           Tract_Had_Redlines +
                           Med_Home_Value_2010_IMP+
                           Ln_Dist_NCAandT_Univ+
                           Ln_NC_CENTRAL_UNIV+
                           Ln_Dist_WSSU+
                           Share_Dem_Vote_2000+
                           Unemp_Rate_2000+
                           Real_GDP_2012+
                           Share_Homeownership_2000+
                           gsmn_math_g3_2013+
                           UNIV_COUNT+
                           Net_Mig_Rate_NWs_1970s+
                           Net_Mig_Rate_NWs_1980s, 
                         data = white_mobility, 
                         method = "lm",
                         metric = "RMSE",
                         trControl = trainControl(method = "cv", number = 10)#,
                         #preProcess = c("zv", "center", "scale")
))

#(white mobility) extract out of sample performance measures
summary(resamples(list(
  Multiple_regression = cv_model1white, 
  MARS = tuned_white
)))$statistics$RMSE %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))

#hispanic mobility

hisp_mobility <- Fina_Dataset_toUse %>%
  select(-c(kir_top20_pooled_pooled_p25,
            kir_top20_black_pooled_p25,
            kir_top20_white_pooled_p25)) %>%
  drop_na()


set.seed(123)
(tuned_hisp <- train(kir_top20_hisp_pooled_p25 ~ hhinc_mean2000+
                        mean_commutetime2000+
                        Share_frac_coll_plus2000+
                        Share_foreign_share2010+
                        popdensity2000+
                        Share_Poor_2000+
                        Share_Black_2000+
                        Share_SingleParent_2000+
                        job_density_2013+
                        Share_Commutes_Less15mins_2010+
                        Share_Commuters_Over1Hr_2010+
                        Share_Commuters_PubTrans_2010+
                        Share_Commuters_PersVeh_2010+
                        Share_HH_withPub_Assist_Inc_2010+
                        Share_Industry_AgForFishMin_2010+
                        Share_Industry_Construction_2010+
                        Share_Industry_Manufact_2010+
                        Share_Industry_RetailTrade_2010+
                        Share_Industry_Information_2010+
                        Share_Industry_Finance_2010+
                        Share_Indus_EdcSrvHelSocAs_2010+
                        Share_Industry_WhlesleTrade_2010+
                        Share_Industry_TrnsWareUtil_2010+
                        Share_Industry_PrfSciManWst_2010+
                        Share_Indus_ArtEntrRecAcc_2010+
                        Share_Indus_PubAdmin_2010+
                        Share_Vacant_Units_2010+
                        Gini_Index_2010+
                        Male_life_expect_2010+
                        Share_Heavy_Drinker_2010+
                        Total_Lynchings_1983to1941+
                        Entropy_Score_2000+
                        `Fe,ale_obesity_prevalence_2009`+
                        Share_Asian_2000_IMP+
                        rent_burden_2010_IMP+
                        Tract_Had_Redlines +
                        Med_Home_Value_2010_IMP+
                        Ln_Dist_NCAandT_Univ+
                        Ln_NC_CENTRAL_UNIV+
                        Ln_Dist_WSSU+
                        Share_Dem_Vote_2000+
                        Unemp_Rate_2000+
                        Real_GDP_2012+
                        Share_Homeownership_2000+
                        gsmn_math_g3_2013+
                        UNIV_COUNT+
                        Net_Mig_Rate_NWs_1970s+
                        Net_Mig_Rate_NWs_1980s,
                      data = hisp_mobility,
                      method = "earth",
                      metric = "RMSE",
                      trControl = trainControl(method = "cv",
                                               number = 10,
                                               savePredictions = "all"),
                      tuneGrid = hyper_grid1
))


tuned_hisp$results %>%
  filter(nprune == tuned_hisp$bestTune$nprune, degree == tuned_hisp$bestTune$degree)


set.seed(123)
(cv_model1hisp <- train(kir_top20_hisp_pooled_p25 ~ hhinc_mean2000+
                           mean_commutetime2000+
                           Share_frac_coll_plus2000+
                           Share_foreign_share2010+
                           popdensity2000+
                           Share_Poor_2000+
                           Share_Black_2000+
                           Share_SingleParent_2000+
                           job_density_2013+
                           Share_Commutes_Less15mins_2010+
                           Share_Commuters_Over1Hr_2010+
                           Share_Commuters_PubTrans_2010+
                           Share_Commuters_PersVeh_2010+
                           Share_HH_withPub_Assist_Inc_2010+
                           Share_Industry_AgForFishMin_2010+
                           Share_Industry_Construction_2010+
                           Share_Industry_Manufact_2010+
                           Share_Industry_RetailTrade_2010+
                           Share_Industry_Information_2010+
                           Share_Industry_Finance_2010+
                           Share_Indus_EdcSrvHelSocAs_2010+
                           Share_Industry_WhlesleTrade_2010+
                           Share_Industry_TrnsWareUtil_2010+
                           Share_Industry_PrfSciManWst_2010+
                           Share_Indus_ArtEntrRecAcc_2010+
                           Share_Indus_PubAdmin_2010+
                           Share_Vacant_Units_2010+
                           Gini_Index_2010+
                           Male_life_expect_2010+
                           Share_Heavy_Drinker_2010+
                           Total_Lynchings_1983to1941+
                           Entropy_Score_2000+
                           `Fe,ale_obesity_prevalence_2009`+
                           Share_Asian_2000_IMP+
                           rent_burden_2010_IMP+
                           Tract_Had_Redlines +
                           Med_Home_Value_2010_IMP+
                           Ln_Dist_NCAandT_Univ+
                           Ln_NC_CENTRAL_UNIV+
                           Ln_Dist_WSSU+
                           Share_Dem_Vote_2000+
                           Unemp_Rate_2000+
                           Real_GDP_2012+
                           Share_Homeownership_2000+
                           gsmn_math_g3_2013+
                           UNIV_COUNT+
                           Net_Mig_Rate_NWs_1970s+
                           Net_Mig_Rate_NWs_1980s, 
                         data = hisp_mobility, 
                         method = "lm",
                         metric = "RMSE",
                         trControl = trainControl(method = "cv", number = 10)#,
                         #preProcess = c("zv", "center", "scale")
))



#(hisp mobility) extract out of sample performance measures
summary(resamples(list(
  Multiple_regression = cv_model1hisp, 
  MARS = tuned_hisp
)))$statistics$RMSE %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))



# principal component regression
set.seed(123)
cv_model2 <- train(kir_top20_pooled_pooled_p25 ~ hhinc_mean2000+
                     mean_commutetime2000+
                     Share_frac_coll_plus2000+
                     Share_foreign_share2010+
                     popdensity2000+
                     Share_Poor_2000+
                     Share_Black_2000+
                     Share_SingleParent_2000+
                     job_density_2013+
                     Share_Commutes_Less15mins_2010+
                     Share_Commuters_Over1Hr_2010+
                     Share_Commuters_PubTrans_2010+
                     Share_Commuters_PersVeh_2010+
                     Share_HH_withPub_Assist_Inc_2010+
                     Share_Industry_AgForFishMin_2010+
                     Share_Industry_Construction_2010+
                     Share_Industry_Manufact_2010+
                     Share_Industry_RetailTrade_2010+
                     Share_Industry_Information_2010+
                     Share_Industry_Finance_2010+
                     Share_Indus_EdcSrvHelSocAs_2010+
                     Share_Industry_WhlesleTrade_2010+
                     Share_Industry_TrnsWareUtil_2010+
                     Share_Industry_PrfSciManWst_2010+
                     Share_Indus_ArtEntrRecAcc_2010+
                     Share_Indus_PubAdmin_2010+
                     Share_Vacant_Units_2010+
                     Gini_Index_2010+
                     Male_life_expect_2010+
                     Share_Heavy_Drinker_2010+
                     Total_Lynchings_1983to1941+
                     Entropy_Score_2000+
                     `Fe,ale_obesity_prevalence_2009`+
                     Share_Asian_2000_IMP+
                     rent_burden_2010_IMP+
                     Tract_Had_Redlines +
                     Med_Home_Value_2010_IMP+
                     Ln_Dist_NCAandT_Univ+
                     Ln_NC_CENTRAL_UNIV+
                     Ln_Dist_WSSU+
                     Share_Dem_Vote_2000+
                     Unemp_Rate_2000+
                     Real_GDP_2012+
                     Share_Homeownership_2000+
                     gsmn_math_g3_2013+
                     UNIV_COUNT+
                     Net_Mig_Rate_NWs_1970s+
                     Net_Mig_Rate_NWs_1980s, 
  data = em_train, 
  method = "pcr",
  trControl = trainControl(method = "cv", number = 10),
  metric = "RMSE",
  preProcess = c("zv", "center", "scale"),
  tuneLength = 20
)

# partial least squares regression
set.seed(123)
cv_model3 <- train(kir_top20_pooled_pooled_p25 ~ hhinc_mean2000+
                     mean_commutetime2000+
                     Share_frac_coll_plus2000+
                     Share_foreign_share2010+
                     popdensity2000+
                     Share_Poor_2000+
                     Share_Black_2000+
                     Share_SingleParent_2000+
                     job_density_2013+
                     Share_Commutes_Less15mins_2010+
                     Share_Commuters_Over1Hr_2010+
                     Share_Commuters_PubTrans_2010+
                     Share_Commuters_PersVeh_2010+
                     Share_HH_withPub_Assist_Inc_2010+
                     Share_Industry_AgForFishMin_2010+
                     Share_Industry_Construction_2010+
                     Share_Industry_Manufact_2010+
                     Share_Industry_RetailTrade_2010+
                     Share_Industry_Information_2010+
                     Share_Industry_Finance_2010+
                     Share_Indus_EdcSrvHelSocAs_2010+
                     Share_Industry_WhlesleTrade_2010+
                     Share_Industry_TrnsWareUtil_2010+
                     Share_Industry_PrfSciManWst_2010+
                     Share_Indus_ArtEntrRecAcc_2010+
                     Share_Indus_PubAdmin_2010+
                     Share_Vacant_Units_2010+
                     Gini_Index_2010+
                     Male_life_expect_2010+
                     Share_Heavy_Drinker_2010+
                     Total_Lynchings_1983to1941+
                     Entropy_Score_2000+
                     `Fe,ale_obesity_prevalence_2009`+
                     Share_Asian_2000_IMP+
                     rent_burden_2010_IMP+
                     Tract_Had_Redlines +
                     Med_Home_Value_2010_IMP+
                     Ln_Dist_NCAandT_Univ+
                     Ln_NC_CENTRAL_UNIV+
                     Ln_Dist_WSSU+
                     Share_Dem_Vote_2000+
                     Unemp_Rate_2000+
                     Real_GDP_2012+
                     Share_Homeownership_2000+
                     gsmn_math_g3_2013+
                     UNIV_COUNT+
                     Net_Mig_Rate_NWs_1970s+
                     Net_Mig_Rate_NWs_1980s, 
  data = em_train,  
  method = "pls",
  trControl = trainControl(method = "cv", number = 10),
  metric = "RMSE",
  preProcess = c("zv", "center", "scale"),
  tuneLength = 20
)

plot(cv_model3)
cv_model3$bestTune
summary(cv_model3)

summary(cv_model3) %>% .$coefficients %>% head(10)

plot(varImp(cv_model3), 14, main = "PLS")

cv_model3$coefnames


# extract out of sample performance measures
summary(resamples(list(
  Multiple_regression = cv_model1, 
  MARS = tuned_mars,
  PLS = cv_model3,
  PCR = cv_model2
)))$statistics$RMSE %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))

cv_model2$preProcess



p1 <- vip(tuned_mars, num_features = 48, geom = "point", value = "gcv") + ggtitle("GCV")

p2 <- vip(tuned_mars, num_features = 48, geom = "point", value = "rss") + ggtitle("RSS")

gridExtra::grid.arrange(p1, p2, ncol = 2)


p1 <- partial(tuned_mars, pred.var = "Share_frac_coll_plus2000", grid.resolution = 10) %>%
  autoplot()
p2 <- partial(tuned_mars, pred.var = "Share_Vacant_Units_2010", grid.resolution = 10) %>%
  autoplot()
p3 <- partial(tuned_mars, pred.var = c("Share_frac_coll_plus2000", "Share_Vacant_Units_2010"), grid.resolution = 10) %>% 
  plotPartial(levelplot = FALSE, zlab = "yhat", drape = TRUE, colorkey = TRUE, screen = list(z = -20, x = -60))

gridExtra::grid.arrange(p1, p2, p3, ncol = 3)



summary(resamples(list(
  Multiple_regression = cv_model1, 
  MARS = tuned_mars
)))$statistics$RMSE %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))

##Gaussian Processes
library(plgp)
n <- 8
X <- matrix(seq(0,2*pi, length=n), ncol = 1)
y <- sin(X)
D <- distance(X)
eps <- sqrt(.Machine$double.eps)
Sigma <- exp(-D) + diag(eps, ncol(D))


XX <- matrix(seq(-0.5, 2*pi+0.5, length=100), ncol = 1)
DXX <- distance(XX)
SXX <- exp(-DXX) + diag(eps, ncol(DXX))
DX <- distance(XX, X)
SX <- exp(-DX)

Si <- solve(Sigma)
mup <- SX %*% Si %*% y
Sigmap <- SXX - SX %*% Si %*% t(SX)


YY <- rmvnorm(100, mup, Sigmap)

q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))

matplot(XX, t(YY),
        type = "l",
        col = "gray",
        lty = 1,
        xlab = "x",
        ylab = "y")
points(X, y, pch=20, cex=2)
lines(XX, sin(XX), col="blue")
lines(XX, mup, lwd=2)
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)

##2D
nx <- 20
x<- seq(0, 2, length=nx)
X <- expand.grid(x, x)
D <- distance(X)
Sigma <- exp(-D) + diag(eps, nrow(X))
Y <- rmvnorm(2, sigma = Sigma)
par(mfrow=c(1,2))
persp(x, x, matrix(Y[1,], ncol = nx), theta = -30, phi = 30, xlab = "x1",
      ylab = "x2", zlab = "y")
persp(x, x, matrix(Y[2,], ncol = nx), theta = -30, phi = 30, xlab = "x1",
      ylab = "x2", zlab = "y")

library(lhs)
X <- randomLHS(40, 2)
