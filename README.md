# SNP를 이용한 품종분류

# 목차

1. 문제

2. snp_information
   - 염색체
   - 센티모르간
   - 위치
   - 유전자 다형성


---

# 문제와 배경지식

유전체 염기서열에서 획득한 유전체 변이 정보인 Single Nucleotide Polymorphism 정보는 특정 개체 및 특정 품종에 따라 다른 변이 양상을 나타낼 수 있기 때문에 동일개체를 확인하거나,
동일 품종을 구분하는데 활용이 가능합니다. 특히 많은 SNP 정보를 통해 분류하는 것보다, 보다 더 적은 SNP 정보로 높은 분류 성능을 내는 것이 중요합니다.

# snp_information

![Example instance](https://github.com/bloodmage1/Snp_approach/blob/main/snp_infp.png)

주어진 정보는 다음과 같습니다. 

## 염색체(Chromosome)

다른 품종은 종종 염색체의 특정 위치에 있는 대립유전자의 빈도로 구별할 수 있다. 전체 게놈에 걸친 일련의 snp를 사용하여 품종별 유전적 변이를 찾을 수 있다.

```
train_x = train_x.drop(["father","mother","gender"],axis = 1)
test_x = test_x.drop(["father","mother","gender"],axis = 1)


for j in range(len(train_x)):
    train_x["SNP_15"][j] = "1" + str(train_x["SNP_15"][j])

for j in range(len(test_x)):
    test_x["SNP_15"][j] = "1" + str(test_x["SNP_15"][j])
    
for j in range(len(train_x)):
    train_x["SNP_02"][j] = "2" + str(train_x["SNP_02"][j])

for j in range(len(test_x)):
    test_x["SNP_02"][j] = "2" + str(test_x["SNP_02"][j])
    
for j in range(len(train_x)):
    train_x["SNP_01"][j] = "3" + str(train_x["SNP_01"][j])
    train_x["SNP_03"][j] = "3" + str(train_x["SNP_03"][j])
    train_x["SNP_10"][j] = "3" + str(train_x["SNP_10"][j])
    train_x["SNP_12"][j] = "3" + str(train_x["SNP_12"][j])
    train_x["SNP_13"][j] = "3" + str(train_x["SNP_13"][j])
    train_x["SNP_14"][j] = "3" + str(train_x["SNP_14"][j])

for j in range(len(test_x)):
    test_x["SNP_01"][j] = "3" + str(test_x["SNP_01"][j])
    test_x["SNP_03"][j] = "3" + str(test_x["SNP_03"][j])
    test_x["SNP_10"][j] = "3" + str(test_x["SNP_10"][j])
    test_x["SNP_12"][j] = "3" + str(test_x["SNP_12"][j])
    test_x["SNP_13"][j] = "3" + str(test_x["SNP_13"][j])
    test_x["SNP_14"][j] = "3" + str(test_x["SNP_14"][j])
    
for j in range(len(train_x)):
    train_x["SNP_04"][j] = "4" + str(train_x["SNP_04"][j])
    train_x["SNP_05"][j] = "4" + str(train_x["SNP_05"][j])
    train_x["SNP_06"][j] = "4" + str(train_x["SNP_06"][j])
    train_x["SNP_07"][j] = "4" + str(train_x["SNP_07"][j])
    train_x["SNP_08"][j] = "4" + str(train_x["SNP_08"][j])
    train_x["SNP_09"][j] = "4" + str(train_x["SNP_09"][j])

for j in range(len(test_x)):
    test_x["SNP_04"][j] = "4" + str(test_x["SNP_04"][j])
    test_x["SNP_05"][j] = "4" + str(test_x["SNP_05"][j])
    test_x["SNP_06"][j] = "4" + str(test_x["SNP_06"][j])
    test_x["SNP_07"][j] = "4" + str(test_x["SNP_07"][j])
    test_x["SNP_08"][j] = "4" + str(test_x["SNP_08"][j])
    test_x["SNP_09"][j] = "4" + str(test_x["SNP_09"][j])

for j in range(len(train_x)):
    train_x["SNP_11"][j] = "5" + str(train_x["SNP_11"][j])

for j in range(len(test_x)):
    test_x["SNP_11"][j] = "5" + str(test_x["SNP_11"][j])

class_le = preprocessing.LabelEncoder()
snp_le = preprocessing.LabelEncoder()
snp_col = [f'SNP_{str(x).zfill(2)}' for x in range(1,16)]

snp_data = []
for col in snp_col:
    snp_data += list(train_x[col].values)
    
train_y = class_le.fit_transform(train_y)
snp_le.fit(snp_data)

for col in train_x.columns:
    if col in snp_col:
        train_x[col] = snp_le.transform(train_x[col])
        test_x[col] = snp_le.transform(test_x[col])

        
train_x["snp_poly"] = use_train["snp_poly"]
test_x["snp_poly"] = use_test["snp_poly"]
```

## 센티모르간(cm)

센티모르간은 유전자 재조합 빈도의 측정 단위이다. 염색체에서 유전자의 위치와 상대적 거리를 mapping하기 위해 유전적 연결 분석에 사용된다.

```
max_value = snp_info["cm"].max()
min_value = snp_info["cm"].min()

minmax_list = []

for i in snp_info["pos"]:
    minmax_list.append(1 - ((i-min_value) / (max_value-min_value)))
    
minmax_list
```


## 위치(position)

많은 단백질은 4차 구조를 가지며 이는 올리고머 분자와 결합하는 여러 폴리펩타이드 사슬로 구성된다. 이러한 단백질의 각 폴리펩타이드 사슬을 서브 유닛이라고 한다. 예를 들어 헤모글로빈은 2개의 알파 서브 유닛과 2개의 베타 서브 유닛으로 구성된다. 각각 4개의 체인에는 헴 포켓이 있는 알파 글로빈 접기가 있다.

```
max_value = snp_info["pos"].max()
min_value = snp_info["pos"].min()

minmax_list = []

for i in snp_info["pos"]:
    minmax_list.append(1 - ((i-min_value) / (max_value-min_value)))
    
minmax_list
```

## 유전적 다형성

일부 다형성은 유전자의 기능을 직접 변경한다. 단백질을 암호화하는 유전자에서 다형성이 발생하면 단백질의 아미노산 서열이 변경되어 기능이 변경될 수 있다. 또한 특정 단백질이 생성되는 양을 변경할 수 있다. 

그러나 다형성의 영향은 상황에 따라 달라, 개인에게 미치는 영향이 다를 수 있다.


```
class_le = preprocessing.LabelEncoder()
snp_le = preprocessing.LabelEncoder()
snp_col = [f'SNP_{str(x).zfill(2)}' for x in range(1,16)]

snp_data = []
for col in snp_col:
    snp_data += list(use_train[col].values)
    
train_y = class_le.fit_transform(train_y)
snp_le.fit(snp_data)

use_train = use_train.drop(["father","mother","gender"],axis = 1)
use_test = use_test.drop(["father","mother","gender"],axis = 1)

use_train["snp_poly"] = 0
use_test["snp_poly"] = 0

count_list = []

for j in range(len(use_train)):
    count = 0
    for col in use_train.columns[1:-1]:
        a =use_train[col].str.split(" ")[j][0]
        b =use_train[col].str.split(" ")[j][1]
        
        if a != b:
            count += 1
    
    count_list.append(count/15)

use_train["snp_poly"] = count_list
    
count_list = []
    
for j in range(len(use_test)):
    count = 0
    for col in use_test.columns[1:-1]:
        a =use_test[col].str.split(" ")[j][0]
        b =use_test[col].str.split(" ")[j][1]
        
        if a != b:
            count += 1
    
    count_list.append(count/15)    


use_test["snp_poly"] = count_list
```