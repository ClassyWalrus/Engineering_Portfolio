@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 15676)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 16212)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 1552)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 16332)

del /F cleanup-ansys-DESKTOP-HOGNH48-16332.bat
