@echo off
copy /Y %1 %2
set /a "ECGkit = %errorlevel%"
del %2
set /a "ECGkit = %ECGkit% | %errorlevel%"
exit %ECGkit%
