@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 6168)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 14732)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 9128)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 11868)

del /F cleanup-ansys-DESKTOP-HOGNH48-11868.bat
