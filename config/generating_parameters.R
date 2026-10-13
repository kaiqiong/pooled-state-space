## Frozen aggregate parameters used to generate the final manuscript simulations.
## These are simulation design constants, not results of the current application run.
## Positive: original FINAL_fitted_models_final.rds, fixed observation variance 400.
## Boundary: original PRIMARY_boundary.rds, fixed observation variance zero.
## Hexadecimal floating literals preserve the original IEEE doubles exactly;
## comments show 17-significant-digit decimal values for readability.
## No participant-level information is stored here.
pss_generating_parameters <- list(
  positive = list(
    beta = c(
      Intercept = 0x1.8d1f18df76b2bp+0, # 1.5512557550016008
      Age = 0x1.6ca8f43bb307p-5, # 0.0445141573123643
      Sex = 0x1.77798e4f46deep-4, # 0.091668659112137368
      BMI = 0x1.22a961b0f3752p-6, # 0.017740579034401911
      Prediabetic = 0x1.89347fdec1827p-4, # 0.095997333042845953
      Diabetic = 0x1.57ed533a98846p-2 # 0.33586626095905492
    ),
    a1 = 0x1.ab6a73f3f744p+5, # 53.426978975280235
    P1 = 0x1.dd62c3fe96d5ap+10, # 1909.5432125542634
    sigma_eta2 = 0x1.3452b6432bcf7p+5, # 38.540386700415119
    sigma_eps2 = 0x1.9p+8 # 400
  ),
  boundary = list(
    beta = c(
      Intercept = 0x1.a417bd9b925ffp+0, # 1.6409872536481769
      Age = 0x1.04422acdc76c6p-5, # 0.031769832230776937
      Sex = 0x1.69d33d81baecap-4, # 0.088336220036292407
      BMI = 0x1.2b105094ebe8ep-5, # 0.036506802927189871
      Prediabetic = 0x1.a0c849b1c5b0ap-5, # 0.050876754698199703
      Diabetic = 0x1.000f2b500dc9ep-2 # 0.25005786586866374
    ),
    a1 = 0x1.c44d7444bb5bp+5, # 56.537819420783649
    P1 = 0x1.b2aec7f638c36p+10, # 1738.730954699929
    sigma_eta2 = 0x1.0e61385a6b549p+5, # 33.797470766449415
    sigma_eps2 = 0x0p+0 # 0
  )
)
