@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 4576)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 10720)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 18876)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 22520)

del /F cleanup-ansys-DESKTOP-HOGNH48-22520.bat
