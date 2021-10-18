
const int camPin = 3;
const int led1Pin = 9;
const int led2Pin = 10;
const int led3Pin = 11;

byte cmd_buffer[3];

int T = 100;
int t_pls = 20;

void setup() {
  pinMode(camPin, OUTPUT);
  pinMode(led1Pin, OUTPUT);
  pinMode(led2Pin, OUTPUT);
  pinMode(led3Pin, OUTPUT);

  Serial.begin(9600);
  Serial.setTimeout(100);
  digitalWrite(led1Pin, HIGH);
}

void loop() {

  if (Serial.available()>0) {
  
    int n = Serial.readBytes(cmd_buffer, 3);
    char cmd = cmd_buffer[0];
    short val = cmd_buffer[1]*256 + cmd_buffer[2]*1;
    if (cmd=='T') {T = val;}
    if (cmd=='L') {
      digitalWrite(led1Pin, val==1);
      digitalWrite(led2Pin, val==2);
      digitalWrite(led3Pin, val==3); 
    }
  }

  digitalWrite(camPin, HIGH);
  delay(t_pls);
  digitalWrite(camPin, LOW);
  delay(T - t_pls);
  
}
