#ifndef XmlWriter_H
#define XmlWriter_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class XmlWriter {

public:

    inline void EndStartedTag(){ fprintf(outfile,">\n"); }
    inline void CloseFile(){ fclose(outfile);printf("XML  file Successfully closed\n"); }
    inline void WriteAttribute(const std::string outAttribute){ Write(" " + outAttribute); }

    void OpenFile(const char *filename);
    void OpenElement(const std::string elementTag);
    void StartOpenTag(const std::string elementTag);
    void CloseEmptyElement();
    void CloseElement(); 
    void WriteData(const std::string outString);
    void FinishFile();


    FILE *outfile;
    int indent;
    int openTags;
    std::vector<std::string> tempOpenTag;

    inline void Write(const std::string str){ fprintf(outfile, "%s", str.c_str());}
    void AddTag(const std::string elementTag);
    void RemoveTag();
};

#endif

