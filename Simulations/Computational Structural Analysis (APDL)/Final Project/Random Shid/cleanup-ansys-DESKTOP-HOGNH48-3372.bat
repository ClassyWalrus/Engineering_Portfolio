@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 16636)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 12784)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 2800)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 3372)

del /F cleanup-ansys-DESKTOP-HOGNH48-3372.bat
