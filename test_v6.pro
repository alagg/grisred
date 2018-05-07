

path = '/data/slam/GREGOR-data/GREGOR-Jun16/gris/09jun16/'
dir='./level0/'
files=path+dir+'09jun16.100'
map_out = './'+'09jun16.100'

dum=rfits_im(files[0],1,dd,hdr)
nrhdr=n_elements(hdr)


;get git revision
gitri=routine_info('gris_v6',/source)
gitfi=file_info(file_dirname(gitri.path)+'/grisred.version')
if gitfi.exists then begin
  gitrev=strarr(2)
  openr,unit,/get_lun,gitfi.name & readf,unit,gitrev & free_lun,unit
endif else gitrev=['n/a','n/a n/a']
print,'GIT-revision: ',gitrev

                                ;add git revision 
      igv=0
      hdr2  = hdr
      hdrstr=string(reform(byte(strjoin(hdr,/single)),80,36*nrhdr))

; ; ; ; ; ; Test making sure that the header satisfies the 80 characters width
;       hdrstr=string(hdrstr,format='(A80)')
;       hdrstr=strtrim(hdrstr,1)


numb    = indgen(10)
alphabet=string(bindgen(1,26)+(byte('a'))[0])
strtmp = ''
for j=0,strlen(gitrev[0])-1 do begin
tmp  = strmid(gitrev[0],j,1)
for i=0,n_elements(numb)-1 do begin
mkstr=where(strmatch(tmp,strtrim(numb[i],1)) eq 1,nmk)
if nmk gt 0 then strtmp +=strtrim(numb[i],1)
endfor ; for i
for i=0,n_elements(alphabet)-1 do begin
mkstr=where(strmatch(tmp,strtrim(alphabet[i],1)) eq 1,nmk)
if nmk gt 0 then strtmp +=strtrim(alphabet[i],1)
endfor ; for i
endfor ; for j
print,strtmp
stop

; ; ; ; WORKED! Test making the systime shorter just keeping the importand info
      tim = BIN_DATE(systime(0))
      tim = string(tim[0],f='(I04)')+'-'+string(tim[1],f='(I02)')+$
            '-'+string(tim[2],f='(I02)')
      sxaddpar,hdrstr,'DATERED',string(tim),$
      ' date of data reduction (yyyy-mm-dd)', before='FILENAME'

; ; ; ; ; Test making sure that the query satisfies the width requirements and 
; ; ; ; ; eliminating any posible character that my cause a problem
      gitrev_tmp=(strsplit(gitrev[0],/extract))[0]
      sxaddpar,hdrstr,'GITREV',gitrev_tmp,$
               ' grisred git revision',before='FILENAME'
      gitrev_tmp=(strsplit(gitrev[1],/extract))[1]  ;Cutting no necessary info
       gitrev_tmp=strmid(gitrev_tmp,strpos(gitrev_tmp,'com')+3)    ; ''
       gitrev_tmp=STRJOIN(STRSPLIT(gitrev_tmp, /EXTRACT,'/'), ' ')
       gitrev_tmp=STRJOIN(STRSPLIT(gitrev_tmp,/regex, /EXTRACT,'.git'))
       gitrev_tmp=string(gitrev_tmp,format='(A18)')
      sxaddpar,hdrstr,'GITREP',gitrev_tmp, $
                ' grisred git repository',before='FILENAME'


; ; ; ; ; ; Test trying to write the GIT query on the comment section 
; ;       gitrev_tmp=(strsplit(gitrev[1],/extract))[1]  ;Cutting no necessary info
;       sxaddpar,hdrstr,'GITREP','git repository',' "'+gitrev_tmp+'"', before='FILENAME',Format='a18'
;       gitrev_tmp=string(gitrev[0])
;       sxaddpar,hdrstr,'GITREV','git revision',' "'+gitrev_tmp+'"',before='FILENAME',Format='a18'


;;;; Test trying to add the git version as a comment
; gitrev_tmp=string(gitrev[0])
; sxaddpar,hdrstr,'COMMENT ', 'grisred git revision'
; sxaddpar,hdrstr,'COMMENT ', string(gitrev_tmp,/PRINT)
; gitrev_tmp=string((strsplit(gitrev[1],/extract))[1])
; sxaddpar,hdrstr,'COMMENT ', 'grisred git repository'
; sxaddpar,hdrstr,'COMMENT ', string(gitrev_tmp,/PRINT)

WRITEFITS, map_out+'cc', data, hdrstr;,append=1

a = readfits(map_out+'cc',hdr2)

end
