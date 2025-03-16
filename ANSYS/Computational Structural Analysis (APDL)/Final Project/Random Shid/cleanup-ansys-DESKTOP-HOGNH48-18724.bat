@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 14840)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 1192)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 808)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 18724)

del /F cleanup-ansys-DESKTOP-HOGNH48-18724.bat
