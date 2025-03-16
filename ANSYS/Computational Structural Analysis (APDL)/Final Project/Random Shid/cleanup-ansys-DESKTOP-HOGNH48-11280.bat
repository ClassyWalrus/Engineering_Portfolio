@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 19032)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 1472)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 7600)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 11280)

del /F cleanup-ansys-DESKTOP-HOGNH48-11280.bat
