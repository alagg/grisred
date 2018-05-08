pro patch_hdr_v6,map_out



print,'--------------------------------------->>>>>>>>>>>>'


;get git revision
gitri=routine_info('gris_v6',/source)
gitfi=file_info(file_dirname(gitri.path)+'/grisred.version')
if gitfi.exists then begin
gitrev=strarr(2)
openr,nlun,/get_lun,gitfi.name & readf,nlun,gitrev & free_lun,nlun
endif else gitrev=['n/a','n/a n/a']
print,'GIT-revision: ',gitrev


data = readfits(map_out+'cc',hdrstr)
tim = BIN_DATE(systime(0))
tim = string(tim[0],f='(I04)')+'-'+string(tim[1],f='(I02)')+$
'-'+string(tim[2],f='(I02)')
sxaddpar,hdrstr,'DATERED',string(tim),$
' date of data reduction (yyyy-mm-dd)', before='FILENAME'

gitrev_tmp=string(gitrev[0])
sxaddpar,hdrstr,'GITREV',gitrev_tmp,$
' grisred git revision',before='FILENAME'

gitrev_tmp=(strsplit(gitrev[1],/extract))[1]  
sxaddpar,hdrstr,'GITREP',gitrev_tmp, $
' grisred git repository',before='FILENAME'

;;remove blank spaces
blank_str=string(bytarr(80)+32B)
zblank=where(hdrstr eq blank_str,nblank)
if nblank gt 0 then for j=0,nblank-1 do hdrstr=remove_blankline(hdrstr)


WRITEFITS, map_out+'cc', data, hdrstr



print,'<<<<<<<<<<<---------------------------------------'

end
