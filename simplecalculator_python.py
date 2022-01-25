print('Calculadora Simples')
l = 1
while l == 1:
  print('Escolha uma opção (1,2,3,4): ')
  print('1. Soma \n 2. Subtração \n 3. multiplicação \n 4. Divisão')

  op = input()
  op = int(op)
  
  for i in list(range(1,5,1)):
    if op == i:
      l = 0
    else:
      l = 1
    if l == 0:
      break
  print('Opção inválida')

print('Digite o primeiro número: ')
n1 = input()
n1 = int(n1)

print('Digite o segundo número: ')
n2 = input()
n2=int(n2)

if op == 1:
  result = n1 + n2
  print('O resultado da soma é: ',result)

if op == 2:
  result = n1 - n2
  print('O resultado da subtração é: ',result)

if op == 3:
  result = n1*n2
  print('O resultado da multiplicação é: ',result)

if op == 4:
  result = n1/n2
  print('O resultado da divisão é: ',result)