@echo off
if "%1"=="" goto error
if not exist "%1" goto error
fpp /ansi /nologo /S. /S.\include\2d.v9 /DBL_NEED_RAND /DBL_NO_FORT_FLUSH /DBL_SPACEDIM=2 /DBL_USE_DOUBLE /DBL_LANG_FORT "%1" | perl strip72 -c > "%1"or
goto end

:error
echo usage: fpp_strip72 {infile}

:end

