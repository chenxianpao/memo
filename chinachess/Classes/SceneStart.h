#ifndef  _SceneStart_H_
#define _SceneStart_H_

#include "cocos2d.h"
USING_NS_CC;

class SceneStart : public CCLayer
{
public:
    static CCScene* scene();

    bool init();

    CREATE_FUNC(SceneStart);

    bool ccTouchBegan(CCTouch* pTouch, CCEvent* pEvent);
    void ccTouchEnded(CCTouch* pTouch, CCEvent* pEvent);

    CCSprite* _red;
    CCSprite* _black;

    void update(float);

    //����Ƿ�ѡ���˺�ɫ����
     bool _selected;
};

#endif

