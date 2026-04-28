@echo off
cd /d "C:\Users\DELL\Desktop\Aerospace\BTech\BTech Y3\C++ CFD\CFD practice\Two Dimensions"
cl.exe /Zi /EHsc /std:c++17 /nologo /Fe:program.exe main.cpp geometry.cpp weno5.cpp rk_ssp.cpp
if %errorlevel% == 0 (
    echo.
    echo Build successful. Running program...
    program.exe
) else (
    echo.
    echo Build failed.
)