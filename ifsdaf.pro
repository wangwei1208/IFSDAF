function cost_fun, g
  common v_pub1,ind_v
  common v_pub2,dep_v
  l=total((ind_v ## g-dep_v)^2)
  return, l
end

pro IFSDAF

;****************************改进的FSDAF算法*************************
compile_opt idl2
envi,/restore_base_save_files
envi_batch_init
;************************************************输入参数***************************************************
;*********************************************************************************************************
fbase_fp="C:\Users\15845\Desktop\NewProject\landsat11.26.dat"                                ;基准landsat
cbase_fp="C:\Users\15845\Desktop\NewProject\modis331.500.dat"                                ;基准modis
cpre_fp="C:\Users\15845\Desktop\NewProject\modis347.500.dat"                                 ;预测modis
fpre_fp="C:\Users\15845\Desktop\NewProject\modis347.500_ifsdaf.dat"                          ;预测结果
em="C:\Users\15845\Desktop\NewProject\SVD_Endmembers.csv"                                    ;端元数据
fc_ratio=20                                                                                  ;modis和landsat像素比例
half_win_size=1                                                                              ;移动窗口半径，值为modis像素，移动窗口大小要大于端元数
ab_threshold=0.05                                                                            ;端元丰度阈值，默认0.05
;**************************************************************************************************************
;******************************************************打开数据**************************************************
envi_open_file,fbase_fp,r_fid=fb_fid
if fb_fid eq -1 then begin
  envi_batch_exit
  return
endif
envi_open_file,cbase_fp,r_fid=cb_fid
if cb_fid eq -1 then begin
  envi_batch_exit
  return
endif
envi_open_file,cpre_fp,r_fid=cp_fid=
if cp_fid eq -1 then begin
  envi_batch_exit
  return
endif
;文件查询
envi_file_query,fb_fid,dims=fb_dims,nl=fb_lines,ns=fb_samples,nb=fb_bands,data_type=fb_dt
envi_file_query,cb_fid,dims=cb_dims,nl=cb_lines,ns=cb_samples,nb=cb_bands,data_type=cb_dt
envi_file_query,cp_fid,dims=cp_dims,nl=cp_lines,ns=cp_samples,nb=cp_bands,data_type=cp_dt
;获取空间参考
fb_mapinfo =envi_get_map_info(fid=fb_fid)
cb_mapinfo =envi_get_map_info(fid=cb_fid)
routine_dir=file_dirname(routine_filepath('IFSDAF'))+'\'  ;获取程序目录

;step 1  时间预测，此过程只能预测土地类型未发生变化的情况
;1 打开端元光谱文件
;*****************************************************************************
em_spec=read_csv(em,header=em_name)
em_samples=n_elements(em_name) 
em_lines  =n_elements(em_spec.(0));em_lines等于波段数
temp=fltarr(em_samples,em_lines)
for i=0,em_samples-1 do temp[i,*]=float(em_spec.(i))
em_spec=temporary(temp);em_spec存储各个波段光谱值

;2 完全最小二乘解混
;*******************************************************************
;移动窗口大小
win_size =2*half_win_size+1
pixelcnts=win_size*win_size
if pixelcnts le (em_samples) then begin ;移动窗口的平方大于端元数
  print,'移动窗口太小'
  envi_batch_exit
  return
endif

;fcls
cd,routine_dir
fbaseabd_fp=file_dirname(fpre_fp)+'\'+file_basename(fbase_fp)+'_abundance.tif'
cmdstr='abundancecaculatemodule.exe '+fbase_fp+' '+em+' '+fbaseabd_fp
spawn,cmdstr,/hide
envi_open_file,fbaseabd_fp,r_fid=fabd_fid
if fabd_fid eq -1 then begin
  envi_batch_exit
  print,'光谱解混失败'
  return
endif

;3 丰度聚集，求取modis像素的丰度
;*********************************************************
;求解光谱角，判断相似端元

lut=uintarr(em_samples,em_samples);存放光谱相似度的索引值，其值大小一次递减
for i=0,em_samples-1 do begin
  cossam=fltarr(em_samples)
  for j=0,em_samples-1 do begin
    if j eq i then continue ;若i=j 进入下一次循环
    cossam[j]=total(em_spec[i,*]*em_spec[j,*])/(sqrt(total(em_spec[i,*]^2))*sqrt(total(em_spec[j,*]^2)))
  endfor
  lut[*,i]=reverse(sort(cossam)) ;lun存放数据值由大到小的索引值
  ;print,"reverse",lut[*,i] 
endfor

;丰度聚集
fabd_img=fltarr(fb_samples,fb_lines,em_samples);存放landsat像素内的丰度值
for ib=0,em_samples-1 do begin
  fabd_img[*,*,ib]=envi_get_data(fid=fabd_fid,pos=ib,dims=fb_dims)
  ;处理异常端元丰度
  index=where(fabd_img lt 0.0,count)
  if count gt 0 then begin
    fabd_img[index,ib]=0.0
  endif
endfor
envi_file_mng,id=fabd_fid,/remove,/delete
cabd_img=fltarr(cb_samples,cb_lines,em_samples)   ;cabd_img存放 modis丰度值
npixels=float(fc_ratio)*fc_ratio ;npixels为一个MODIS像素中有多少个Landsat

;小于阈值的丰度合并到最相似的端元上
for ci=0,cb_lines-1 do begin
  for cj=0,cb_samples-1 do begin
    ;逐MODIS像素遍历
    win_data=fabd_img[(cj*fc_ratio):((cj+1)*fc_ratio-1),(ci*fc_ratio):((ci+1)*fc_ratio-1),*] ;win_data存储的一个MODIS像素内所有Landsat像素的丰度值
    cabd=fltarr(em_samples);cabd存储的临时粗像素的丰度，也就是从landsat聚集到modis中
    for ib=0,em_samples-1 do begin
      cabd[ib]=mean(win_data[*,*,ib])
    endfor
    index=where((cabd gt 0 ),cnts) ;此处未考虑一个modis像素内多个端元丰度小于0.05
    abd_min=min(cabd[index],min_ind) ;min_ind返回最小值的索引，abd_min的值为最小丰度值
    while abd_min lt ab_threshold do begin
      min_ind=index[min_ind] ;将cabd[index]中的索引转为cabd中的索引
      abdtemp=cabd[min_ind]
      cabd[min_ind]=0.0
      ;将小于0.05的端元赋予最相似的端元   
      for sc=0,em_samples-1 do begin
        if cabd[lut[sc,min_ind]] ne 0 then begin ;第一列为最相似的端元
          cabd[lut[sc,min_ind]]+=abdtemp
          break
        endif
      endfor       
      index=where((cabd gt 0),cnts)
      abd_min=min(cabd[index],min_ind);最小的值小于0.05已经加到最相似的端元了，若第二小也小于0.05,急需下一次循环，若最小的值不满足while循环，退出
    endwhile
    cabd_img[cj,ci,*]=cabd
  endfor
endfor 

;4 约束最小二乘解混，求取modis像素内每一类的变化
;*************************************************************
change_f =make_array(fb_samples,fb_lines,cb_bands,type=cb_dt)  
cb_img   =make_array(cb_samples,cb_lines,cb_bands,type=cb_dt)
cp_img   =make_array(cp_samples,cp_lines,cp_bands,type=cp_dt)
fb_img   =make_array(fb_samples,fb_lines,fb_bands,type=fb_dt)
for nb=0,cb_bands-1 do cb_img[*,*,nb]=envi_get_data(fid=cb_fid,dims=cb_dims,pos=nb)
for nb=0,cp_bands-1 do cp_img[*,*,nb]=envi_get_data(fid=cp_fid,dims=cp_dims,pos=nb)
for nb=0,fb_bands-1 do fb_img[*,*,nb]=envi_get_data(fid=fb_fid,dims=fb_dims,pos=nb)
envi_file_mng,id=cb_fid,/remove
envi_file_mng,id=cp_fid,/remove
envi_file_mng,id=fb_fid,/remove
change_c=cp_img-cb_img
cf_img=make_array(cb_samples,cb_lines,cb_bands,type=cb_dt);存储的是基准日期一个modis像素内所有landsat的平均值
he_index=fltarr(fb_samples,fb_lines,fb_bands)

;设置端元的变化阈值 
min_allow=fltarr(em_samples,cb_bands)
max_allow=fltarr(em_samples,cb_bands)
for ib=0,cb_bands-1,1 do begin
  min_allow[*,ib]=min(change_c[*,*,ib])-stddev(change_c[*,*,ib])
  max_allow[*,ib]=max(change_c[*,*,ib])+stddev(change_c[*,*,ib])
endfor
;约束最小二乘参数设置
common v_pub1
common v_pub2
gbnd    =[0,100]     
nobj    = 0           
lcomp   = 'cost_fun'  
for ci=0,cb_lines-1 do begin
  for cj=0,cb_samples-1 do begin
    ai=max([0,ci-half_win_size])       ; 移动窗口直径
    bi=min([cb_lines-1,ci+half_win_size])
    aj=max([0,cj-half_win_size])
    bj=min([cb_samples-1,cj+half_win_size])
    fai=ci*fc_ratio
    fbi=(ci+1)*fc_ratio-1
    faj=cj*fc_ratio
    fbj=(cj+1)*fc_ratio-1 ;landsat窗口大小
    c_win_pixels=(bi-ai+1)*(bj-aj+1)  
;    print,c_win_pixels  
    if (c_win_pixels gt em_samples) then begin   ;窗口大小大于端元数，继续执行
      fabd_temp=fabd_img[fai:fbi,faj:fbj,*]        ;窗口内landsat像素存在fabd_temp    
;      cabd_temp=cabd_img[ci,cj,*]
      ind_v   = transpose(reform(cabd_img[aj:bj,ai:bi,*],c_win_pixels,em_samples))
      for nb=0,cb_bands-1 do begin
         dep_v = double(reform(change_c[aj:bj,ai:bi,nb],c_win_pixels))
         x     = fltarr(1,em_samples)
         xbnd  = [[min_allow[*,nb]], [max_allow[*,nb]]]
         constrained_min, x, xbnd, gbnd, nobj, lcomp,inform, nstop = 5 ;此处x返回的是每个modis像素内的端元变化范围
         ds_change=fabd_temp*rebin(reform(x,1,1,em_samples),fc_ratio,fc_ratio,em_samples,/sample);最邻近采样为了保证值不变的情况下进行矩阵相乘
         change_f[faj:fbj,fai:fbi,nb]=total(ds_change,3);change_f为每一landsat的变化值
         cf_img[cj,ci,nb]=mean(fb_img[faj:fbj,fai:fbi,nb]);cf_img为一landsat像素聚集为modis像素的值，便于调整传感器差异         
      endfor
    endif    
  endfor
endfor
;5 传感器差异调整，求取系数变化
;*********************************************************************
fp_img_temp=make_array(fb_samples,fb_lines,fb_bands,type=fb_dt)
for ib=0,fb_bands-1  do begin
  x =reform(cb_img[*,*,ib],cb_samples*cb_lines)
  y =reform(cf_img[*,*,ib],cb_samples*cb_lines)
 ; print,"x",x
;  print,"y",y
  coefs=linfit(x,y)
  fp_img_temp[*,*,ib]=fb_img[*,*,ib]+(coefs[1]*change_f[*,*,ib])
endfor
;时间预测值进行聚合
cp_img_temp=fltarr(cb_samples,cb_lines,cb_bands)
for ci=0,cb_lines-1 do begin
  for cj=0,cb_samples-1 do begin
    for ib=0,fb_bands-1 do begin
      cp_img_temp[cj,ci,ib]=mean(fp_img_temp[cj:((cj+1)*fc_ratio-1),ci:((ci+1)*fc_ratio-1),ib])
    endfor    
  endfor  
endfor
;  
;    
       
;6 时间预测，求取最终预测值
;***********************************************************************
openw,lun,fpre_fp,/get_lun
writeu,lun,fp_img_temp
free_lun,lun
envi_setup_head,fname=fpre_fp,nb=fb_bands, ns=fb_samples,nl=fb_lines, interleave=0,data_type=fb_dt,map_info=fb_mapinfo,/write

;step 2  空间预测 


;7 TPS薄板函数插值，求取土地覆盖发生变化信息
;************************************************************************
;计算landsat 行列索引数组
row_ind=intarr(fb_samples,fb_lines)
col_ind=intarr(fb_samples,fb_lines)
for i=0,fb_samples-1,1 do begin
  col_ind[i,*]=i
endfor
for i=0,fb_lines-1,1 do begin
  row_ind[*,i]=i
endfor
;计算modis行列索引数组
row_c=fltarr(cb_samples,cb_lines)
col_c=fltarr(cb_samples,cb_lines)
for ci=0,cb_lines-1 do begin
  for cj=0,cb_samples-1 do begin
    row_c[cj,ci]=mean(row_ind[cj:(cj+1)*fc_ratio-1,ci:(ci+1)*fc_ratio-1])
    col_c[cj,ci]=mean(col_ind[cj:(cj+1)*fc_ratio-1,ci:(ci+1)*fc_ratio-1])
  endfor
endfor

;tps插值
tps=fltarr(fb_samples,fb_lines,fb_bands)
for ib=0,fb_bands-1,1 do begin
  tps[*,*,ib] = min_curve_surf(cp_img[*,*,ib], col_c, row_c,/tps, xpout=col_ind,  ypout=row_ind)
endfor
;fn=dialog_pickfile(title="tps结果")
;openw,lun,fn,/get_lun
;writeu,lun,tps
;free_lun,lun
;envi_setup_head,fname=fn,nb=fb_bands, ns=fb_samples,nl=fb_lines, interleave=0,data_type=fb_dt,map_info=fb_mapinfo,/write

;*******************************************************************************************
;**********************************************************************************************
;计算fb_img中的相似像素的阈值,后面进行计算同质系数以及加权求和要用
similar_th=fltarr(fb_bands)
for iband=0,fb_bands-1,1 do begin
  similar_th[iband]=stddev(fb_img[*,*,iband])*2.0/em_samples
endfor
;同质系数计算
he_index=fltarr(fb_samples,fb_lines,fb_bands)
for fi=0,fb_lines-1 do begin
  for fj=0,fb_samples-1 do begin
    ai=max([0,fi-(half_win_size*fc_ratio)])
    bi=min([fb_lines-1,fi+half_win_size*fc_ratio]) ;搜索窗口，[fj,fi]为中心像素    
    aj=max([0,fj-(half_win_size*fc_ratio)])
    bj=min([fb_samples-1,fj+half_win_size*fc_ratio])
;    print,"lajissssssssssssssssssssssssssss",ai,bi,aj,bj
    temp=fltarr(bj-aj+1,bi-ai+1,fb_bands)
    for ib=0,fb_bands-1 do begin
      temp[*,*,ib]=fb_img[aj:bj,ai:bi,ib]-fb_img[fj,fi,ib]
      index=where(temp[*,*,ib] lt similar_th[ib],nums)
      if nums gt 0 then begin
        he_index[fj,fi,ib]=nums/((bj-aj+1)*(bi-ai+1))
      endif else begin
        he_index[fj,fi,ib]=0.0
      endelse      
    endfor
  endfor 
endfor

;;*********************************************************************************************
;********************************************************************************************
;8 残差分配
;**************************************************************************
predict_change_c=cp_img_temp-cf_img     
real_change_c=cp_img-cb_img
change_r=real_change_c-predict_change_c ;change_r为残差值
change_21_c=fltarr(cb_samples,cb_lines,cb_bands) ;残差分配
change_21=fltarr(fb_samples,fb_lines,fb_bands);landsat每一像素残差分配结果
for ci=0,cb_lines-1 do begin
  for cj=0,cb_samples-1 do begin
    fai=ci*fc_ratio
    fbi=(ci+1)*fc_ratio-1
    faj=cj*fc_ratio
    fbj=(cj+1)*fc_ratio-1 ;modis像素内landat的索引
    fb_nums=float(fc_ratio)*fc_ratio ;一个modis内多少个lansat像素
    for ib=0,cb_bands-1 do begin
      diff_change=change_r[cj,ci,ib]
      w_change_tps=(tps[*,*,ib])[faj:fbj,fai:fbi]-(fp_img_temp[*,*,ib])[faj:fbj,fai:fbi];文中Eh0
      if (diff_change le 0) then begin ;diff-change 小于0，则w_change_tps也要小于 0
        ind_noc=where(w_change_tps gt 0, num_noc)
        if (num_noc gt 0) then begin
          w_change_tps[ind_noc]=0
        endif
      endif else begin
        ind_noc=where(w_change_tps lt 0, num_noc)
        if (num_noc gt 0) then begin
          w_change_tps[ind_noc]=0
        endif
      endelse
      w_change_tps=abs(w_change_tps)
      w_unform=fltarr(fb_nums)     ;对于landsat像素级别
      w_unform[*]=abs(diff_change)
      w_change=w_change_tps*he_index[faj:fbj,fai:fbi,ib]+w_unform*(1.0-he_index[faj:fbj,fai:fbi,ib])+0.000001  ;combine these two weights
      w_change=w_change/(mean(w_change)) ;归一化
      ;去除异常值
      ind_extrem=where(w_change gt 10, num_extrem)
      if (num_extrem gt 0) then begin
        w_change[ind_extrem]=mean(w_change)
      endif
      w_change=w_change/(mean(w_change)) ;w_change为加权值
      change_21[faj:fbj,fai:fbi,ib]=w_change*diff_change
    endfor
  endfor
endfor
;将change_21加到fp_img_temp上
;允许t2变化范围
min_allow=fltarr(fb_bands)
max_allow=fltarr(fb_bands)
for ib=0,fb_bands-1 do begin
  min_allow[ib]=min([min(cp_img[*,*,ib]),min(fp_img_temp[*,*,ib])])
  max_allow[ib]=max([max(cp_img[*,*,ib]),max(fp_img_temp[*,*,ib])])
endfor

fp_img_redis=fltarr(fb_samples,fb_lines,fb_bands)
fp_img_redis=fP_img_temp+change_21
for ib=0,fb_bands-1 do begin
  temp=fP_img_redis[*,*,ib]
  index_min=where(temp lt min_allow[ib],min_nums)
  if min_nums gt 0 then begin
    temp[index_min]=min_allow[ib]
  endif
  index_max=where(temp gt max_allow[ib],max_nums)
  if max_nums gt 0 then begin
    temp[index_max]=max_allow[ib]
  endif
  fp_img_redis[*,*,ib]=temp
  change_21[*,*,ib]=fp_img_redis[*,*,ib]-fp_img_temp[*,*,ib] ;残差值
endfor
;fpre_fp1="C:\Users\15845\Desktop\ifsdaf.dat"
;openw,lun,fpre_fp1,/get_lun
;writeu,lun,fp_img_redis
;free_lun,lun
;envi_setup_head,fname=fpre_fp1,nb=fb_bands, ns=fb_samples,nl=fb_lines, interleave=0,data_type=fb_dt,map_info=fb_mapinfo,/write


;9 加权函数进行最终预测  
;d_d_all为窗口内到中心像素的距离
fp_img=fltarr(fb_samples,fb_lines,fb_bands);最终预测结果
w=half_win_size*fc_ratio;landsat像素
d_d_all=((w-indgen(w*2+1)#(intarr(1,w*2+1)+1))^2+(w-(intarr(w*2+1)+1)#indgen(1,w*2+1))^2)^0.5 ;d_d_all为完整窗口中心像素到周围其他像素的距离
d_d_all=reform(d_d_all,(w*2+1)*(w*2+1))
;第一次写的
;for fi=0,fb_lines-1 do begin
;  for fj=0,fb_samples-1 do begin
;      fai=max([0,fi-w])
;      fbi=min([fb_lines-1,fi+w])
;      faj=max([0,fj-w])
;      fbj=min([fb_samples-1,fj+w]) ;窗口范围   
;      i_tar=fi-fai
;      j_tar=fj-faj ;目标像素位置
;      temp_similar=fltarr((fbj-faj+1),(fbi-fai+1)) ;存储同一窗口内与目标像素的差值
;      similar=fltarr((fbj-faj+1),(fbi-fai+1),fb_bands) ;窗口内像素与目标像素的度量
;      postion_similar=intarr((fbj-faj+1),(fbi-fai+1),fb_bands);存放相似像素的位置
;      ;寻找相似像素
;      for ib=0,fb_bands-1 do begin
;        win_img=fb_img[faj:fbj,fai:fbi,ib]
;        temp_similar=win_img-win_img[j_tar,i_tar]
;        similar[*,*,ib]=similar[*,*,ib]+abs(temp_similar)/win_img[j_tar,i_tar];值越小，代表相似度越高
;        index=where(abs(temp_similar)lt similar_th[ib],similar_nums ) ;寻找相似像        
;        postion_similar[index,ib]=1 ;将满足条件的像素赋值为1
;      endfor
;      ;如果窗口不完整
;      col_wind=indgen(fbi-fai+1)#(intarr(1,fbj-faj+1)+1)*1.0 ;列的值都一样
;      row_wind=(intarr(fbi-fai+1)+1)#indgen(1,fbj-faj+1)*1.0 ;行的值都一样
;      for ib=0,fb_bands-1 do begin
;        similar_index=where(postion_similar[*,*,ib] eq 1 ,nums)
;        if  ((fbi-fai+1)*(fbj-faj+1) lt (w*2.0+1)*(w*2.0+1)) then begin           
;            d_d=((col_wind[similar_index])^2+(row_wind[similar_index])^2)^0.5 ;距离公式
;        endif else begin
;            d_d=d_d_all[similar_index]
;        endelse
;        d_d=(1.0+d_d/w)*(similar[similar_index,ib]+1.0)
;        c_d=1.0/d_d
;        weight=c_d/total(c_d)
;        ;最终预测
;        change_cand=(change_21[faj:fbj,fai:fai,ib])[similar_index]
;        fp_img[j_tar,i_tar,ib]=total(weight*change_cand)+fb_img[j_tar,i_tar,ib]
;      endfor      
;  endfor
;endfor


;第二次写的
for fi=0,fb_lines-1 do begin
  for fj=0,fb_samples-1 do begin
    fai=max([0,fi-w])
    fbi=min([fb_lines-1,fi+w])
    faj=max([0,fj-w])
    fbj=min([fb_samples-1,fj+w]) ;窗口范围
    i_tar=fi-fai
    j_tar=fj-faj ;目标像素位置    
    col_wind=indgen(fbi-fai+1)#(intarr(1,fbj-faj+1)+1)*1.0 ;列的值都一样
    row_wind=(intarr(fbi-fai+1)+1)#indgen(1,fbj-faj+1)*1.0 ;行的值都一样
    ;search similar pixels within window
    similar_cand=fltarr((fbi-fai+1)*(fbj-faj+1)) ;pleace the similarity measure between each pixel and the target pixel
    position_cand=intarr((fbi-fai+1)*(fbj-faj+1))+1  ;place the location of each similar pixel
    for ib=0,fb_bands-1,1 do begin
      cand_band=intarr((fbi-fai+1)*(fbj-faj+1))
      wind_fine=fb_img[faj:fbj,fai:fbi,ib]
      s_s=abs(wind_fine-wind_fine[j_tar,i_tar])
      similar_cand=similar_cand+s_s/(wind_fine[j_tar,i_tar]+0.00000001)
      ind_cand=where(s_s lt similar_th[ib])
      cand_band[ind_cand]=1
      position_cand=position_cand*cand_band
    endfor
    indcand=where(position_cand ne 0,number_cand)  ;select similar pixel initially
    order_dis=sort(similar_cand[indcand])
    ;number_cand=min([number_cand0,num_similar_pixel])
    ind_same_class=indcand[order_dis[0:number_cand-1]]           ; select the n most similar samples
    ;compute weight for each simialr pixel
    ;spatial distance
    if ((fbi-fai+1)*(fbj-faj+1) lt (w*2.0+1)*(w*2.0+1)) then begin   ;not an integrate window
      d_d_cand=((i_tar-col_wind[ind_same_class])^2+(j_tar-row_wind[ind_same_class])^2)^0.5+0.00001
    endif else begin
      d_d_cand=d_d_all[ind_same_class]      ;integrate window
    endelse

    ;normalize these distances
    d_d_cand=(1.0+d_d_cand/w)*(similar_cand[ind_same_class]+1.0)
    c_d=1.0/d_d_cand
    weight=c_d/total(c_d)
    w_pixel=float(fc_ratio)*fc_ratio
    for ib=0,fb_bands-1,1 do begin
      change_cand=(change_21[faj:fbj,fai:fbi,ib])[ind_same_class]
      fp_img[fj,fi,ib]=fb_img[fj,fi,ib]+total(weight*change_cand)
      x=fb_img[fj,fi,ib]+cp_img[ceil(fj/w_pixel),ceil(fi/w_pixel),ib]-cb_img[ceil(fj/w_pixel),ceil(fi/w_pixel),ib]
      if fp_img[fj,fi,ib] lt 0.0 then begin
        another_predict=max([0,x])
        fp_img[fj,fi,ib]=min([1.0,another_predict])        
      endif
      if fp_img[fj,fi,ib] gt 1.0 then begin
        another_predict=min([1.0,x])
        fp_img[fj,fi,ib]=max([0,another_predict])
      endif
    endfor
  endfor
endfor
      

fpre_fp1="C:\Users\15845\Desktop\ifsdaf_result.dat"
openw,lun,fpre_fp1,/get_lun
writeu,lun,fp_img
free_lun,lun
envi_setup_head,fname=fpre_fp1,nb=fb_bands, ns=fb_samples,nl=fb_lines, interleave=0,data_type=fb_dt,map_info=fb_mapinfo,/write  
    



















end