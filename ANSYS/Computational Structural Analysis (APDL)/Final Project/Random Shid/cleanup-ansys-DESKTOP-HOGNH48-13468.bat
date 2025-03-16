@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 14580)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 11456)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 2740)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 13468)

del /F cleanup-ansys-DESKTOP-HOGNH48-13468.bat
