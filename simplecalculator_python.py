print('Calculadora Simples')
l = 1
while l == 1:
  print('Escolha uma op��o (1,2,3,4): ')
  print('1. Soma \n 2. Subtra��o \n 3. multiplica��o \n 4. Divis�o')

  op = input()
  op = int(op)
  
  for i in list(range(1,5,1)):
    if op == i:
      l = 0
    else:
      l = 1
    if l == 0:
      break
  print('Op��o inv�lida')

print('Digite o primeiro n�mero: ')
n1 = input()
n1 = int(n1)

print('Digite o segundo n�mero: ')
n2 = input()
n2=int(n2)

if op == 1:
  result = n1 + n2
  print('O resultado da soma �: ',result)

if op == 2:
  result = n1 - n2
  print('O resultado da subtra��o �: ',result)

if op == 3:
  result = n1*n2
  print('O resultado da multiplica��o �: ',result)

if op == 4:
  result = n1/n2
  print('O resultado da divis�o �: ',result)