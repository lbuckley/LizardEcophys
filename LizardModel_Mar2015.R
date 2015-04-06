  # INCLUDE BIOPHYSICAL MODELS
 16 Nov, 2011

#BIOPHYSICAL MODELING ***************************************
# Calculates operative environmental temperature for Sceloporus undulatus in the US (from Buckley 2008, AmNat)

# read climate data


#10' DATA
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Lab\\SceloporusModel\\")

climdata=read.csv('ClimateData_2160x684_noNA.csv');

#read surface temperature data
tsurfdata=read.csv("Tsurface_NARRall_nohead.csv");

lat= climdata[,3] ; #latitude in degrees
lon= climdata[,2] ; #longitude in degrees

Tm= tsurfdata[,2:13]; #mean of an average day of each month
Tr= tsurfdata[,14:25]; #diurnal temperature range of an average day of each month
Wind= climdata[,4]; #mean wind speed m/s
Albedo= climdata[,5:8]; #Albedo percentage for Jan / Apr / Jul / Oct
Elev= climdata[,33]; #mean elevation (m)
TSmax=climdata[,34:45];  #max soil T, mean 14:00hr temperature for five days in the middle of each month (K), from LDAS
TSr=climdata[,46:57];  #range between 02:00 and 14:00hr temperature for five days in the middle of each month (K), from LDAS

Topt=33 #optimal temperature
#__________________________________________________
# Biophysical models from Campbell & Norman 1998
# constants
sigma=5.67*10^-8 # stefan-boltzman constant
c_p=29.3 # specific heat of air, J/mol degrees K or C

# absorptivity
alpha_S=0.9  # solar absorptivity (Gates 1980, Table 11.4)
alpha_L=0.965 # thermal absoptivity, Bartlett & Gates 1967 
epsilon_s=0.965 # surface emisivity (p163), Bartlett & Gates 1967

F_d=0.8  # diffuse view angle, Bartlett & Gates 1967
F_r=0.5  # reflected solar radiation
F_a=0.5  # atmospheric radiation
F_g=0.5  # ground thermal radation

tau=0.65 # atmospheric transmisivity
S_p0=1360 # extraterrestrial flux density (p159)

rd=180/pi  # factor to convert radians into degrees

#--------------------------------------------
#LIZARD DATA

lizdat=read.csv("POPmeans_16Nov2011.csv") 

svl= lizdat[, 8] #mean of max male and female SVL (Ord & Blumstein 2002, Herrel 2002)
mass= 3.55*10^-5*(svl)^3.00 #Tinkle and Ballinger 1972 (g)

m=lizdat[,13]  #mean lizdat[,18], spec: lizdat[,9]
mu=lizdat[,12]  #mean lizdat[,17], spec: lizdat[,8]
a=lizdat[,11] #a mean

# scel undu
VTmin=  lizdat[, 15] 
VTmax=  lizdat[, 16]
VTmean=lizdat[, 2] #PBT

CTmin=  lizdat[, 3]
CTmax=  lizdat[, 4]

eis= lizdat[,14]*0.76*0.5  #assume catch 50# of insects

#--------------------------------------------
# data structures

K= matrix(NA, nrow(Tm),2) #carrying capacity array
Kmat= array(NA, dim=c(nrow(Tm),nrow(lizdat),2)) #carrying capacity matrix for all species
LMat=matrix(NA, nrow(Tm),2)

#Thermoreg then no thermoreg in last dimension
VMat= array(NA, dim=c(nrow(Tm),12,2) ) # velocity matrix
EMat=array(NA, dim=c(nrow(Tm),12,2) ) # digestion efficiency matrix
CMat=array(NA, dim=c(nrow(Tm),12,2) ) # Consumption matrix
ForageMat=array(NA, dim=c(nrow(Tm),12,2) ) # foraging matrix

Tf = matrix(NA, nrow(Tm),2) # time foraging

TeMat=array(NA, dim=c(nrow(Tm),12,2) ) 
TeMatS=array(NA, dim=c(nrow(Tm),12,2) ) 
ForageM=array(NA, dim=c(nrow(Tm),nrow(lizdat),2) )  # foraging matrix

#--------------------------------------------
# FIXED PARAMETERS
L = 1000 # ei=72.3*0.76*0.10 #energetic input per insect (J) #assume catch 10%
T=24*30*60*60 #total number seconds each month

#--------------------------------------------
# Loop through species
#for(Spec in 1:2){ #cell counter

#SPECIFY SPECIES
Spec=1

# ENERGETICS	#regressions derived from Gillooly et al. 2001
                	# assume max MR= 3* RMR from Nagy 2005		               

#ew=lizdat[Spec,17]*1.5 #data for nj and south carolina
ew= exp(-10.0+0.51*log(mass[Spec])+0.115*VTmean[Spec]) *1.5;  #Angilletta PBZ 2001, MR j/s
#ew= lizdat[Spec,7]/(60*60)*1.5;  #POP DATA
#ew=lizdat[Spec,17]*1.5 #data for nj and south carolina

ep= 3/1.5* ew; #calculate maximal energy in joules/sec, factor from bennett 1982

# CALCULATE VELOCITIES OVER TEMPERATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
								# assume forage at 70% of max velocity (Irschick and Losos 1998)                     
vmax=10^(0.044+0.20*log10(mass[Spec]))*0.7;                 # van damme 2001, across lizards
                                
ei=eis[Spec];

#------------------------------------
# Calculate Operative Environmental Temperatures 
#for (cellk in 10000:11000){
for (cellk in 1:nrow(Tm)){ # cell counter
     for (monthk in 1:12){ # month counter
         ForageHrs_T=0; ForageHrs_NT=0
         vtot_T=0; vtot_NT=0
         etot_T=0; etot_NT=0
         Cm_T=0; Cm_NT=0
         dhourk=0; # daylight hour counter

	for (dayk in 1:30){    # day counter
             J= (monthk-1)*30+dayk # Julian calendar day             
             RevAng = 0.21631 + 2 * atan (0.967 * tan (0.0086 * (-186 + J))) # Revolution angle in radians
             DecAng = asin (0.39795 * cos (RevAng))                         # Declination angle in radians           
  
             f=(279.575+0.9856*J)/rd  # f in radians
             ET= (-104.7*sin (f)+596.2*sin (2*f)+4.3*sin (3*f)-12.7*sin (4*f)-429.3*cos (f)-2.0*cos (2*f)+19.3*cos (3*f))/3600   # (11.4) Equation of time
             LC= 1/15* (15 - lon[cellk]%%15) # longitude correction, 1/15h for each degree e of standard meridian
             t_0  =12-LC-ET # solar noon 

             latr = lat[cellk]/rd # latitude in radians
             Daylength = 24 - (24 / pi) * acos ((sin (6 * pi / 180) + sin (latr) * sin (DecAng)) / (cos (latr) * cos (DecAng)))
    
         for (hourk in 1:24) {   # hour counter              
             cpsi_rad= sin (DecAng)*sin (latr) + cos (DecAng)*cos (latr)*cos (pi/12*(hourk-t_0)) # zenith angle in radians
             psi=acos (cpsi_rad) # (11.1) zenith angle in radians
             if (psi>pi/2) psi=pi/2
             
             #__________________________________________________ 
                # Calculate the diurnal temperature trend, Campbell and Norman 1998
                W=pi/12
                gamma= 0.44 - 0.46* sin (0.9 + W * hourk)+ 0.11 * sin (0.9 + 2 * W * hourk)   # (2.2) diurnal temperature function
                Tx= Tm[cellk, monthk] + Tr[cellk, monthk] / 2  # maximum daily temperature
                Tn= Tm[cellk, monthk] - Tr[cellk, monthk] / 2 # minimum daily temperature
                Ta  = Tx+Tn * (1-gamma) 
                 
                # can add temperature increase (e.g., 3C) to simulate climate change          
                #__________________________________________________
                # Calculate radiation
                # view angles, parameterize for animal suspended above ground (p181), on ground- adjust F_e, F_r, and F_g
                h=svl[Spec]/1000 # length of cylinder in m
                theta = psi # angle between solar beam and a normal to the plane in radians, = psi for horizontal surfaces

                # F_p=(cos (theta)+(4*h*sin (theta))/(pi*d))/(2+4*h/d)  # beam view angle, Fig 11.6
                A=0.121*mass[Spec]^0.688   # total lizard area, roughgarden 1981 from Norris (1965) and Porter and James (1979)
                A_p= (-1.1756810^-4*psi^2-9.2594*10^-2*psi+26.2409)*A/100      # projected area
                F_p=A_p/A
                
                # radiation
                albk= floor (monthk/3.01)+1 # compute season for albedo measurement
                rho_S= Albedo[cellk, albk]/100 # (Table 12.2)

                p_a=101.3* exp (-Elev[cellk]/8200)  # atmospheric pressure
                m_a=p_a/(101.3*cos (psi))  # (11.12) optical air mass
                if(psi>80*pi/180) m_a=5.66
                
                # Surface Soiltemp (8.5), degrees K
                # W (8.7) angular frequency for diurnal fluctuation 
                TSa= TSmax[cellk, monthk]-0.5*TSr[cellk, monthk] # Average soil temp (K)
                As= 0.5*TSr[cellk, monthk] # Amplitude of diurnal temperature flux
                Ts= TSa + As*sin (W*(hourk-8))  # t=8 is phase adjustment, see p25
                
                # Flux densities
                epsilon_ac= 9.2*10^-6*(Ta+273)^2 # (10.11) clear sky emissivity
                L_a=epsilon_ac*sigma*(Ta+273)^4  # (10.7) long wave flux densities from atmosphere 
                L_g=epsilon_s*sigma*(Ts)^4  # (10.7) long wave flux densities from ground

                S_p=S_p0*tau^m_a # (11.11) direct irradience 
                S_d=0.3*(1-tau^m_a)* S_p0*cos (psi)   # (11.13) diffuse radiation
                S_t=S_p*cos (psi)+S_d # solar irradience 
                S_r= rho_S*S_t # (11.10) reflected radiation

                R_abs= alpha_S*(F_p*S_p+ F_d*S_d + F_r*S_r)+alpha_L*(F_a*L_a+F_g*L_g) # (11.14) Absorbed radiation
                #__________________________________________________
                
                # conductance
                dim=svl[Spec]/1000 # characteristic dimension in meters (Table 9.5)
                g_r= 4*sigma*(Ta+273)^3/c_p # (12.7) radiative conductance
                g_Ha=1.4*0.135*sqrt(Wind[cellk]/dim) # boundary conductance, factor of 1.4 to account for increased convection (Mitchell 1976)
                
                #__________________________________________________
                # operative environmental temperature
                Te=Ta+(R_abs-epsilon_s*sigma*(Ta+273)^4)/(c_p*(g_r+g_Ha))                       
			#TeMat[cellk, monthk]=Te;
         
                # calculate in shade, set absorbed radiation to 0# of ambient
                TeS=Ta+(0.6*R_abs-epsilon_s*sigma*(Ta+273)^4)/(c_p*(g_r+g_Ha))
         		#TeMatS[cellk, monthk]=TeS;
			#Assume 60% less radiation in shade
                 
day=FALSE
if(hourk>t_0-round(Daylength/2) && hourk<t_0+round(Daylength/2)) day=TRUE #daylight window in hours

if(day==TRUE){ #daylight window in hours, examine every hour
                    
#THERMOREGULATION
#Assume thermoregulate to Topt if possible
if(TeS<Topt && Te>Topt) To=Topt
if(!(TeS<Topt && Te>Topt)) if(TeS<CTmax[Spec]) To=runif(1,min=TeS, max=min(Te, CTmax[Spec]))
if(TeS>=CTmax[Spec]) To=TeS

if(To> VTmin[Spec] && To< VTmax[Spec]){ #function within operative environmental temperatures
                    ForageHrs_T = ForageHrs_T + 1; #foraging
                            if(To> 28.4) vfact= 95 +(40.3-28.4)/5*(To-28.4)
				    if(To<=28.4) vfact= 80 +(28.4-23.03)/15*(To-23.03)
                            vtot_T=vtot_T+vmax*(vfact/100) } #end check operative temperature

#Calculate digestion
                    arcde=85.34-0.50*To+0.000074*To^3; #digestive efficiency from Grant and Porter 1992
                    etot_T= sin(pi/180*arcde)^2  +etot_T; 
                    
      # Maximum consumption from Angilletta 2001 in Joules/day
            if(To<=20)    Cmax=94;
            if (To>20 && To<=30)  Cmax= 94 +(270-94)/10*(To-20);
            if (To>30 && To<=33)  Cmax= 270 +(511-270)/3*(To-30);
            if (To>33 && To<=36) Cmax= 511-(511-421)/3*(To-33);
            if (To>36)  Cmax= 421;
	#Calculate consumption
		Cm_T= Cmax*mass[Spec]+ Cm_T;

#NO THERMOREGULATION
if(TeS<CTmax[Spec]) To= runif(1,TeS, min(Te, CTmax[Spec]))
if(TeS>=CTmax[Spec]) To=TeS

if(To> VTmin[Spec] && To< VTmax[Spec] ){ #function within operative environmental temperatures
                    ForageHrs_NT = ForageHrs_NT + 1; #foraging
                            if(To> 28.4) vfact= 95 +(40.3-28.4)/5*(To-28.4)
				    if(To<=28.4) vfact= 80 +(28.4-23.03)/15*(To-23.03)
                            vtot_NT=vtot_NT+vmax*(vfact/100) 
} #end check operative temperature

#Calculate digestion
                    arcde=85.34-0.50*To+0.000074*To^3; #digestive efficiency from Grant and Porter 1992
                    etot_NT= sin(pi/180*arcde)^2  +etot_NT; 
                    
      # Maximum consumption from Angilletta 2001 in Joules/day
            if(To<=20)    Cmax=94;
            if (To>20 && To<=30)  Cmax= 94 +(270-94)/10*(To-20);
            if (To>30 && To<=33)  Cmax= 270 +(511-270)/3*(To-30);
            if (To>33 && To<=36) Cmax= 511-(511-421)/3*(To-33);
            if (To>36)  Cmax= 421;
	#Calculate consumption
		Cm_NT= Cmax*mass[Spec]+ Cm_NT;

			
} #END CHECK DAY
	  
if(day==FALSE){ #night time
                    R_abs=0; #night time no radiation
                    To=Ta+(-epsilon_s*sigma*(Ta+273)^4)/(c_p*(g_r+g_Ha)); 
               
#Calculate digestion
                    arcde=85.34-0.50*To+0.000074*To^3; #digestive efficiency from Grant and Porter 1992
                    etot_NT= sin(pi/180*arcde)^2  +etot_NT; 
			  etot_T= sin(pi/180*arcde)^2  +etot_T;
                    
      # Maximum consumption from Angilletta 2001 in Joules/day
            if(To<=20)    Cmax=94;
            if (To>20 && To<=30)  Cmax= 94 +(270-94)/10*(To-20);
            if (To>30 && To<=33)  Cmax= 270 +(511-270)/3*(To-30);
            if (To>33 && To<=36) Cmax= 511-(511-421)/3*(To-33);
            if (To>36)  Cmax= 421;
	#Calculate consumption
		Cm_T= Cmax*mass[Spec]+ Cm_T; 
		Cm_NT= Cmax*mass[Spec]+ Cm_NT;

		 } #END CHECK NIGHT 
                    
} #end hours
} #end day

ForageMat[cellk, monthk,]= c(ForageHrs_T, ForageHrs_NT)
         if(vtot_T>0) VMat[cellk, monthk,1]=vtot_T/ForageHrs_T;
	   if(vtot_NT>0) VMat[cellk, monthk,2]=vtot_NT/ForageHrs_NT;
         EMat[cellk, monthk,]= c(etot_T/(24*30),etot_NT/(24*30))
         CMat[cellk, monthk,]= c(Cm_T/(24*30),Cm_NT/(24*30))
} # end month loop
print(cellk)
} # end cell loop

#________________________________________________
 Tf_T=(rowSums(ForageMat[,,1])*60*60)/(12*30) #daily foraging time in seconds
 Vmean_T= rowMeans(VMat[,,1], na.rm=TRUE)
 Emean_T= rowMeans(EMat[,,1], na.rm=TRUE)
 Cmean_T= rowMeans(CMat[,,1], na.rm=TRUE)

 Tf_NT=(rowSums(ForageMat[,,2])*60*60)/(12*30) #daily foraging time in seconds
 Vmean_NT= rowMeans(VMat[,,2], na.rm=TRUE)
 Emean_NT= rowMeans(EMat[,,2], na.rm=TRUE)
 Cmean_NT= rowMeans(CMat[,,2], na.rm=TRUE)

 #________________________________________________      
T=24*60*60 #total number seconds in 24 hours
for(cellk in 1:length(Tm[,6])){
#for (cellk in 10000:11000){
  
# CALCULATE CLUMPED PARAMETERS    
        # adjust ei for digestive efficiency
	#Thermoreg
	et_T=ei*Emean_T[cellk] #adjust ei for digestive efficiency
      v_T=Vmean_T[cellk]
      b_T= m[Spec] * Tf_T[cellk]   
      nu_T= mu[Spec] + m[Spec] * (T-Tf_T[cellk]) * ew  
	#No thermoreg
	et_NT=ei*Emean_NT[cellk] #adjust ei for digestive efficiency
	v_NT=Vmean_NT[cellk]
      b_NT= m[Spec] * Tf_NT[cellk]   
      nu_NT= mu[Spec] + m[Spec] * (T-Tf_NT[cellk]) * ew 
#______________________________________________________
# Calculate Max Consumption
  # Calculate cutoff radius, rs (from Roughgarden 1997, p313) 
rs_T= (-ep+ew+sqrt((ep-ew)^2 +a[Spec]*v_T*(et_T^2)))/(a[Spec]*et_T)
l_T= (-Cmean_T[cellk]*v_T-Tf_T[cellk]*v_T*ew)/(a[Spec]* rs_T*(Cmean_T[cellk]*rs_T-Tf_T[cellk]*v_T*et_T+rs_T*Tf_T[cellk]*ep))
rs_NT= (-ep+ew+sqrt((ep-ew)^2 +a[Spec]*v_NT*(et_NT^2)))/(a[Spec]*et_NT)
l_NT= (-Cmean_NT[cellk]*v_NT-Tf_NT[cellk]*v_NT*ew)/(a[Spec]* rs_NT*(Cmean_NT[cellk]*rs_NT-Tf_NT[cellk]*v_NT*et_NT+rs_NT*Tf_NT[cellk]*ep))

if(!is.na(l_T))if(l_T>1)l_T=1
if(!is.na(l_T))if (l_T<0)l_T=1
if(!is.na(l_NT))if(l_NT>1)l_NT=1
if(!is.na(l_NT))if (l_NT<0)l_NT=1

LMat[cellk,]=c(l_T, l_NT)

 # CALCULATE K   
 # With thermoreg 
  	if (!is.na(b_T) && !is.na(nu_T))if (b_T>0 && nu_T>0 ){  
         q1= b_T*ep+nu_T
         q2= b_T*ew+nu_T
         deta= b_T^2*et_T^2*(a[Spec]*l_T)^2*v_T^2
         detb= 4*a[Spec]*l_T*v_T*(nu_T+b_T*ep)*(nu_T+b_T*ew)
         det= sqrt(deta-detb)
      	q3= b_T*et_T*a[Spec]*l_T*v_T+ det
         K[cellk,1]=(L*q3)/(2*q2*v_T)
         if(!is.na(K[cellk,1]))if(K[cellk,1]<1) K[cellk,1]=0
	} #end 
         
	if (!is.na(b_T) && !is.na(nu_T)) if (b_T<=0 || nu_T<=0) K[cellk,1]=0

 # No thermoreg
  	if (!is.na(b_NT) && !is.na(nu_NT))if (b_NT>0 && nu_NT>0 ){  
         q1= b_NT*ep+nu_NT
         q2= b_NT*ew+nu_NT
         deta= b_NT^2*et_NT^2*(a[Spec]*l_NT)^2*v_NT^2
         detb= 4*a[Spec]*l_NT*v_NT*(nu_NT+b_NT*ep)*(nu_NT+b_NT*ew)
         det= sqrt(deta-detb)
      	q3= b_NT*et_NT*a[Spec]*l_NT*v_NT+ det
         K[cellk,2]=(L*q3)/(2*q2*v_NT)
         if(!is.na(K[cellk,2]))if(K[cellk,2]<1) K[cellk,2]=0
	} #end 
         
	if (!is.na(b_NT) && !is.na(nu_NT)) if (b_NT<=0 || nu_NT<=0) K[cellk,2]=0

} #end loop cells
 
Kmat[,Spec,]=K[,]
ForageM[,Spec,]=cbind(Tf_T, Tf_NT)
#}  # end species loop
 
  #######################################################   

# EXPORT CARRYING CAPACITIES
GridID= climdata[,1]
M=cbind(GridID, lon, lat, Kmat[,,1],Kmat[,,2])
write.csv(M, "ScelPred_30Nov2011.csv") 

#------------------------------------------
#Make map

library(lattice) #for plotting
library(sp)
library(spatstat)
library(latticeExtra)
library(fields)
library(maps)
library(maptools)
library(classInt)
library(colorRamps)

dat<-read.csv("ScelPred_30Nov2011_clean.csv")
#Make spatial points data frame
xy.sp= cbind(dat$lon, dat$lat)
xy.cc= coordinates(xy.sp)
bbox(xy.sp)

#Make Spatial pixels data frame
#grd1 <- SpatialPixelsDataFrame(points=xy.sp, data = dat[,4:5], tolerance=0.1, proj4string=CRS("+proj=longlat +proj=lcc"))
grd1 <- SpatialPixelsDataFrame(points=xy.sp, data = dat[,4:5], tolerance=0.1)

#Make sp.layout
wrld.m= map("world", plot=FALSE, fill=FALSE)
wrld.p= map2SpatialLines(wrld.m, proj4string=CRS("+proj=longlat +proj=lcc"))  #turn map into spatial lines
wrld<- list("sp.lines", wrld.p, col="darkgray")

#K_T
#Define quantiles
quants<-round(quantile(grd1$K_T, probs = seq(0, 1, length.out=49), na.rm=TRUE), digits=3)
quants<-unique(quants)

#specify colors
#cols<-heat.colors(length(quants)+1, alpha=1)
#cols<-cols[length(cols):1]
cols<-tim.colors(length(quants)+1)

txt1= list("sp.text", c(-105, 40), "K with thermoregulation", cex=1.2)
move.layout<-list(wrld, txt1)

#get quantiles and fix legend
labelat = seq(0,length(cols),by=5) 
labs<- levels(as.factor(round(quants,3)))[labelat] #NEED TO FIX
KT=spplot(grd1, "K_T", at=quants, col.regions=cols, sp.layout=move.layout, colorkey=list(at=1:length(quants),labels = list(at=labelat, labels=labs )))

txt1= list("sp.text", c(-105, 40), "K without thermoregulation", cex=1.2)
move.layout<-list(wrld, txt1)
KNT=spplot(grd1, "K_NT", at=quants, col.regions=cols, sp.layout=move.layout, colorkey=list(at=1:length(quants),labels = list(at=labelat, labels=labs )))

#PLOT 2 TOGETHER 
setwd("\\\\Bioark.bio.unc.edu\\buckleylauren\\Work\\ScelUnduPhysiology\\Figs\\Figs_30Nov2011\\")
pdf("Kmap_30Nov2011.pdf",height = 5, width = 12)
par(mfrow=c(2,1), tcl=-0.5,mar=c(0.5,0.5,0.5,0.5), cex.lab=1.1, cex.axis=1, lwd=1, las=1, bty="l", pch=21, cex=100);
plot(KT, split=c(1,1,2,1), more=TRUE)
plot(KNT, split=c(2,1,2,1), more=TRUE)
dev.off()

