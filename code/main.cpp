#include <iostream>
#include "HoleFiller.h"
#include <spdlog/sinks/stdout_sinks.h>
//#include "MeshSmooth.hpp"
//#include "MeshSmooth.h"
#include "MeshDoctor.h"
#include "RBF_interpolation.hpp"
#include <SplineInterpolation.hpp>
#include "FindInnerBoundary.h"
using namespace std;
using namespace MeshFix;

void TransferJsonToVertical(std::string fileName);
//TransferJsonToVertical("C:/Users/st/Desktop/ToothRelated/AliSegmentData/data/lower_label.txt");
//1.spdlog ��־��Դ��
//2.OpenMesh �ǽṹ����Դ��
//3.Hermite ��ֵ�Ҽ��ɵõ�HoleFiller������
int main(int argc, char*argv[])
{
	Triangle_mesh inputMesh;
	bool inputFlag = OpenMesh::IO::read_mesh(inputMesh, "../../DataDir/lower.stl");//����������Ŀ¼�ҵ����������ļ���
	if (!inputFlag)
		return -1;
	FindInnerBoundary findBoundary(inputMesh);
	findBoundary.SetFaceLabelMap("../../DataDir/lower_label.txt");
	findBoundary.TraverseFaces();
	return system("pause");
}

//��json�ļ�ת��Ϊ��ֱ��
void TransferJsonToVertical(std::string fileName)
{
	ifstream ifile(fileName);
	//���ļ����뵽ostringstream����buf��
	ostringstream buf;
	char ch;
	while (buf && ifile.get(ch))
	{
		switch (ch)
		{
		case ',':
			continue;
		case ' ':
			ch = '\n';
		default:
			buf.put(ch);;
		}
		
	}
	//������������buf�������ַ���
	//std::string a = buf.str();
	std::ofstream m_OutFile("D:/ttttttt.txt");
	m_OutFile << buf.str();
	m_OutFile.close();
	
}