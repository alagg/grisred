pro cal09jun16

cd,'/data/slam/GREGOR-data/GREGOR-Jun16/gris/09jun16/'
lambda=10830.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fileff=['09jun16.001']
filecal=['09jun16.000']
map='09jun16.002'
gris_v6,map,fileff,filecal,lambda=lambda,/show
gris_cc2fits4d,'./level1/09jun16.002*cc'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fileff=['09jun16.003','09jun16.005']
filecal=['09jun16.008']
map='09jun16.004'
gris_v6,map,fileff,filecal,lambda=lambda
gris_cc2fits4d,'./level1/09jun16.004*cc'

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fileff=['09jun16.005','09jun16.007']
filecal=['09jun16.008']
map='09jun16.006'
gris_v6,map,fileff,filecal,lambda=lambda,/show
gris_cc2fits4d,'./level1/09jun16.006*cc'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fileff=['09jun16.009','09jun16.011']
filecal=['09jun16.008']
map='09jun16.010'
gris_v6,map,fileff,filecal,lambda=lambda,/show
gris_cc2fits4d,'./level1/09jun16.010*cc'

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fileff=['09jun16.012','09jun16.014']
filecal=['09jun16.015']
map='09jun16.013'
gris_v6,map,fileff,filecal,lambda=lambda,/show
gris_cc2fits4d,'./level1/09jun16.013*cc'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fileff=['09jun16.017','09jun16.019']
filecal=['09jun16.016']
map='09jun16.018'
gris_v6,map,fileff,filecal,lambda=lambda,/show
gris_cc2fits4d,'./level1/09jun16.018*cc'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

return
end
