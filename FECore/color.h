#include <windows.h>
#include <stdio.h>

void SetConsoleColour(WORD* Attributes, DWORD Colour)
{
	CONSOLE_SCREEN_BUFFER_INFO Info;
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
	GetConsoleScreenBufferInfo(hStdout, &Info);
	*Attributes = Info.wAttributes;
	SetConsoleTextAttribute(hStdout, Colour);
}

void ResetConsoleColour(WORD Attributes)
{
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), Attributes);
}