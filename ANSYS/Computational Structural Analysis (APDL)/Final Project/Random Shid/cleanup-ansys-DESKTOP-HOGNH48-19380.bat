@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 12232)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 7012)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 10556)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 19380)

del /F cleanup-ansys-DESKTOP-HOGNH48-19380.bat
