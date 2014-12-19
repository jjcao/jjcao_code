@Echo Off

@Echo Find debug

@For /r . %%a In (.) Do @If Exist "%%a\debug" @Echo "%%a\debug"

@Echo Find debug Dir....OK

@For /r . %%a In (.) Do @If Exist "%%a\debug" rd /s /q "%%a\debug"

@Echo Find release

@For /r . %%a In (.) Do @If Exist "%%a\release" @Echo "%%a\release"

@Echo Find release Dir....OK

@For /r . %%a In (.) Do @If Exist "%%a\release" rd /s /q "%%a\release"

@Echo Find ipch

@For /r . %%a In (.) Do @If Exist "%%a\ipch" @Echo "%%a\ipch"

@Echo Find ipch Dir....OK

@For /r . %%a In (.) Do @If Exist "%%a\ipch" rd /s /q "%%a\ipch"


@Echo Clear debug Dir Mission Completed

del *.ncb /s
del *.sdf /s

@echo on