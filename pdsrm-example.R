################################
# This code accompanies Knight & Humphrey's chapter on dyadic data analysis. 

# To cite this code, please use:
#  Knight, A. P., & Humphrey, S. E. (2019). Dyadic data analysis. In S. E. Humphrey and J. M. LeBreton (Eds.), The Handbook of Multilevel Theory, Measurement, and Analysis, pp. 423-447. Washington, DC: American Psychological Association.


# The code conducts a social relations analysis on a dataset using a new method
# for the lme function. 
################################

source("http://apknight.org/pdSRM.R")					# This provides the pdSRM functions needed to estimate the SRM

#################################
# The following models were used to produce the results in Knight & Humphrey
#################################

## Null model for cognitive trust
cog.0 <- 
lme(trust_cog ~ 										# Focal criterion variable, which is a directed dyadic rating
	1,													# Explicitly specifying the intercept term
	random = list(										# Begin to specify the random effects portion as a list
		team_id = pdBlocked(list(						# Blocked structure for team_id, with a list containing two elements
			pdIdent(~1), 								# Group-level intercept
			pdSRM(~-1 									# Individual-level part, -1 means no additional group-level intercept	
				+ a1 + a2 + a3  + a4 					# Individual-level (cont.), this is the actor part
				+ p1 + p2 + p3 + p4)))),				# Individual-level (cont.), this is the partner part
	correlation=corCompSymm(form=~1 | team_id/dyad_id), # Dyad-level part, with a compound symmetric structure
	data=d.sub, na.action=na.omit)							# Give the input data name, remove observations with missing data
	
summary(cog.0)											# Summarize the output, which gives parameter estimates
srm.pct(cog.0)											# Output the variance parameters and variance decomposition

cog.1 <- 

lme(trust_cog ~ 									# Focal criterion variable, which is a directed dyadic rating
1													# Explicitly specifying the intercept term
+ gender_pct_grd + social_skills_x_grd				# Group-level fixed effects
+ act_gender + act_social_skills_grd				# Individual-level actor fixed effects
+ part_gender + part_social_skills_grd				# Individual-level partner fixed effects
+ act_gender*part_gender + absdif_social_skills_grd,# Dyad-level fixed effects
random = list(										# Begin to specify the random effects portion as a list
team_id = pdBlocked(list(						# Blocked structure for team_id, with a list containing two elements
pdIdent(~1), 								# Group-level intercept
pdSRM(~-1 									# Individual-level part, -1 means no additional group-level intercept	
+ a1 + a2 + a3  + a4 					# Individual-level (cont.), this is the actor part
+ p1 + p2 + p3 + p4)))),				# Individual-level (cont.), this is the partner part
correlation=corCompSymm(form=~1 | team_id/dyad_id), # Dyad-level part, with a compound symmetric structure
data=d.sub, na.action=na.omit)							# Give the input data name, remove observations with missing data
	
summary(cog.1)											# Summarize the output, which gives parameter estimates
srm.pct(cog.1)											# Output the variance parameters and variance decomposition


aff.0 <- lme(trust_aff ~ 
			1,
			random = list(
				team_id = pdBlocked(list(
					pdIdent(~1), 
					pdSRM(~-1 + a1 + a2 + a3  + a4 
					+ p1 + p2 + p3 + p4)))),				
			correlation=corCompSymm(form=~1 | team_id/dyad_id), 
		data=d.sub, na.action=na.omit)
summary(aff.0)
srm.pct(aff.0)

aff.1 <- lme(trust_aff ~ 
			gender_pct_grd + social_skills_x_grd
			+ act_gender + act_social_skills_grd
			+ part_gender + part_social_skills_grd
			+ act_gender*part_gender + absdif_social_skills_grd
			,
			random = list(
				team_id = pdBlocked(list(
					pdIdent(~1), 
					pdSRM(~-1 + a1 + a2 + a3  + a4 
					+ p1 + p2 + p3 + p4)))),				
			correlation=corCompSymm(form=~1 | team_id/dyad_id), 
		data=d.sub, na.action=na.omit)
summary(aff.1)
srm.pct(aff.1)

