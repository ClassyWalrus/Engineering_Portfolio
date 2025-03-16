@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 12720)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 13944)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 12340)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 2732)

del /F cleanup-ansys-DESKTOP-HOGNH48-2732.bat
