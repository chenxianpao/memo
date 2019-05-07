#ifndef  _Stone_H_
#define  _Stone_H_

#include "cocos2d.h"
USING_NS_CC;


//������
class Stone : public CCSprite
{
public:

    Stone();

    //���ӵ�����:����ʿ���ࡢ�������ڡ���
    enum TYPE {JIANG,SHI,XIANG,CHE,MA,PAO,BING};

    //��������
    //��һ�����������ӵ�����
    //�ڶ������������ӵ���ɫ
    static Stone* create(int id, bool red);
    
    //��ʼ������
     bool init(int id, bool red);

     //������
    void reset(bool red);

    //����÷�����ӵĳ�ʼλ��
    static struct InitPos
    {
        int _x;
        int _y;
       Stone::TYPE _type;
    }_initPos[16];


     CC_SYNTHESIZE(TYPE, _type, Type)
    /*
    protected:
     Type _type;
public:
    void setType(Type var);
    {
        _type = type;
    }
    Type getType() const
    {
        return _type;
    }
*/

    CC_SYNTHESIZE(int, _x, X)
    CC_SYNTHESIZE(int, _y, Y)
    CC_SYNTHESIZE(int, _id, ID)
    CC_SYNTHESIZE(bool, _dead, Dead)
    CC_SYNTHESIZE(bool, _red, Red)
   /*//���ӵ�λ��(����)
    int x;
    int y;

    int _id;//���ӵ�ID  0~31(һ����32������)

    bool dead;//�ж������Ƿ񱻳���

    bool _red;//�ж����ӵ���ɫ*/
};

#endif
