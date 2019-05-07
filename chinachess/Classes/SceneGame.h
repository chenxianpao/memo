#ifndef  _SceneGame_H_
#define _SceneGame_H_

#include "cocos2d.h"
#include "SimpleAudioEngine.h"
#include "Stone.h"
USING_NS_CC;
using namespace CocosDenshion;


//����һ����
//��¼ÿ��һ�������Ϣ
class Step : public CCObject
{
public:
    int _moveid;//��Ҫ�ƶ������ӵ�id
    int _killid;//ͨ���������λ���жϴ��������Ƿ�������
    int _xFrom; //���ӵ�ǰ��λ�õ�x����
    int _yFrom; //���ӵ�ǰ��λ�õ�y����
    int _xTo;   //�����ƶ����λ�õ�x����
    int _yTo;   //�����ƶ����λ�õ�y����
    
    static Step* create(int moveid, int killid, int xFrom, int yFrom, int xTo, int yTo)
    {
        Step* step = new Step;

        step->_moveid = moveid;
        step->_killid = killid;
        step->_xFrom = xFrom;
        step->_yFrom = yFrom;
        step->_xTo = xTo;
        step->_yTo = yTo;

        step->autorelease();

        return step;
    }
};

class SceneGame : public CCLayer
{
public:
     SceneGame();

    ~SceneGame()
    {
        _steps->release();
    }

    static CCScene* scene(bool red);

    //�Զ���init����
    bool init(bool red);

    //�Զ���create����
    //red�����ж���ҽ�����Ϸʱʱѡ�е����ӵ���ɫ
    //�����ѡ�к���ʱ����ҳֺ���
    //�����ѡ�к���ʱ����ҳֺ���
    static SceneGame*  create(bool red);

    //���̵�ƫ����
    CCPoint _plateOffset;
    
    //���ӵ�ƫ����
    CCPoint _stoneOffset; 
    
    //���ӵ�ֱ��
    float _d;

    //����Ƿ��ߺ���
    bool _redTrun;
    
    //������ҵ������Ǻ�ɫ���Ǻ�ɫ
    bool _redSide;
    
    //���ڴ���ѡ���(��������)
    CCSprite* _selectSprite;

    //����ѡ�е����ӵ�id
    int _selectid;

     //����ÿ���ߵ�����
    CCArray* _steps;

     //�������Ӷ�������(������һ����32������)
    Stone* _s[32];

    //���ڴ���������ʾ��Ϸ���
    CCSprite* sprite;

    //�ж���Ϸ�������ʾ״̬
    bool visible;

    //���ð���ʱ,���ӵ�λ��
    void SetRealPos(Stone* s);

    //�õ�������������ϵ������
    //���������λ���������ⷵ��false
    //ͨ��������������������
    bool getClickPos(CCPoint ptInWin, int &x, int &y);

    //ͨ��������±��ȡ���ӵ�ID
    int getStone(int x, int y);

    //ѡ������
    void setSelectId(int id);

    //�ƶ�����
    //��һ���������ƶ�������
    //�ڶ�����������ɱ��������
    void moveStone(int moveId, int killId, int x, int y);

    //�����̵�����ת���ɴ��ڵ�����
    CCPoint getStonePos(int x, int y);

    //����(xo,yo)��(x,y)֮���������
    //���������Ϊ-1,��ʾ(xo,yo)��(x,y)����һ��ֱ����
    int getStoneCount(int xo, int yo, int x, int y);

    //�¾�
    void New(CCObject*);

    //����
    void Back(CCObject*);

    //��ʼ��Ϸ
    void Start(CCObject*);

    //��ͣ��Ϸ
    void Pause(CCObject*);

    //������Ϸ�Ѷ�
    void Difficulty(CCObject*);

    //���ű�������
     void Voice(CCObject*);

    void moveComplete(CCNode*, void*);

    //�������
    bool canMove(int moveid, int killid, int x, int y);

    //�����������
    bool canMoveJiang(int moveid, int killid, int x, int y);

    //ʿ���������
    bool canMoveShi(int moveid, int x, int y);

    //����������
    bool canMoveXiang(int moveid, int x, int y);

    //�����������
    bool canMoveChe(int moveid, int x, int y);

    //����������
    bool canMoveMa(int moveid, int x, int y);

    //�ڵ��������
    bool canMovePao(int moveid, int killid, int x, int y);

    //�����������
    bool canMoveBing(int moveid, int x, int y);
          
    //ͨ�����ѡ�����ӣ�������
    bool ccTouchBegan(CCTouch* pTouch, CCEvent* pEvent);
};

#endif // SCENEGAME_H
