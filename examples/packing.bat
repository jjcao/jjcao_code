@Echo Off

del *.ncb /s
del *.sdf /s
del *.filters /s
del *.user /s
del *.suo /s


@Echo Find debug

@For /r . %%a In (.) Do @If Exist "%%a\debug" @Echo "%%a\debug"

@Echo Find debug Dir....OK

@For /r . %%a In (.) Do @If Exist "%%a\debug" rd /s /q "%%a\debug"

@Echo Clear debug Dir Mission Completed

@Echo Find release

@For /r . %%a In (.) Do @If Exist "%%a\release" @Echo "%%a\release"

@Echo Find debug Dir....OK

@For /r . %%a In (.) Do @If Exist "%%a\release" rd /s /q "%%a\release"

@Echo Clear debug Dir Mission Completed

@echo on