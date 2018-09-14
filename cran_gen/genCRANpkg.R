path = getwd()
path = substr(path, start=1, stop=nchar(path)-9)
print(path)
dir_data = paste0(path,'/data')
unlink(dir_data, recursive = TRUE)

src_dir = paste0(path,'/R-proj/src')
src_dist = path
file.copy(src_dir, src_dist, recursive=TRUE)

inst_dir = paste0(path,'/R-proj/inst')
inst_dist = path
file.copy(inst_dir, inst_dist, recursive=TRUE)

man_dir = paste0(path,'/R-proj/man')
man_dist = path
file.copy(man_dir, man_dist, recursive=TRUE)

R_dir = paste0(path,'/R-proj/R')
R_dist = path
file.copy(R_dir, R_dist, recursive=TRUE)

tests_dir = paste0(path,'/R-proj/tests')
tests_dist = path
file.copy(tests_dir, tests_dist, recursive=TRUE)

descr_dir = paste0(path,'/R-proj/DESCRIPTION')
descr_dist = path
file.copy(descr_dir, descr_dist, recursive=TRUE)

namesp_dir = paste0(path,'/R-proj/NAMESPACE')
namesp_dist = path
file.copy(namesp_dir, namesp_dist, recursive=TRUE)

volpro_dir = paste0(path,'/R-proj/volesti.Rproj')
volpro_dist = path
file.copy(volpro_dir, volpro_dist, recursive=TRUE)

external_dir = paste0(path,'/external')
external_dist = paste0(path,'/src')
file.copy(external_dir, external_dist, recursive=TRUE)

include_dir = paste0(path,'/include')
include_dist = paste0(path,'/src')
file.copy(include_dir, include_dist, recursive=TRUE)

#circleci = paste0(path,'/.circleci')
#unlink(circleci, recursive = TRUE)

readme_dir = paste0(path,'/cran_gen/README.md')
readme_dist = path
file.copy(readme_dir, readme_dist, recursive=TRUE)

news_dir = paste0(path,'/cran_gen/NEWS.md')
news_dist = path
file.copy(news_dir, news_dist, recursive=TRUE)

#makefile_dir = paste0(path,'/cran_gen/Makefile')
#makefile_dist = paste0(path, '/src/lp_solve')
#file.copy(makefile_dir, makefile_dist, recursive=TRUE)

#makevars_dir = paste0(path,'/cran_gen/Makevars')
#makevars_dist = paste0(path, '/src')
#file.copy(makevars_dir, makevars_dist, recursive=TRUE)

#makevarswin_dir = paste0(path,'/cran_gen/Makevars.win')
#makevarswin_dist = paste0(path, '/src')
#file.copy(makevarswin_dir, makevarswin_dist, recursive=TRUE)

#volume_dir = paste0(path,'/cran_gen/volume.h')
#volume_dist = paste0(path,'/src/include/volume'
#file.copy(volume_dir, volume_dist, recursive=TRUE)

#dir_ext = paste0(path,'/external')
#unlink(dir_ext, recursive = TRUE)

#dir_inc = paste0(path,'/include')
#unlink(dir_inc, recursive = TRUE)

#dir_boost = paste0(path,'/src/external/boost')
#unlink(dir_boost, recursive = TRUE)

#dir_eigen = paste0(path,'/src/external/eigen')
#unlink(dir_eigen, recursive = TRUE)



