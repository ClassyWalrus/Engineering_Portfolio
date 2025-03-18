@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 23484)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 13736)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 22664)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 9240)

del /F cleanup-ansys-DESKTOP-HOGNH48-9240.bat
