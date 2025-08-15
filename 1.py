import sqlite3


connection = sqlite3.connect("university.db")
cursor = connection.cursor()


cursor.execute('''
CREATE TABLE IF NOT EXISTS Students (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    surname TEXT NOT NULL,
    age INTEGER NOT NULL,
    city TEXT NOT NULL
);
''')


cursor.execute('''
CREATE TABLE IF NOT EXISTS Courses (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    time_start TEXT NOT NULL,
    time_end TEXT NOT NULL
);
''')


cursor.execute('''
CREATE TABLE IF NOT EXISTS Student_courses (
    student_id INTEGER NOT NULL,
    course_id INTEGER NOT NULL,
    FOREIGN KEY (student_id) REFERENCES Students (id),
    FOREIGN KEY (course_id) REFERENCES Courses (id)
);
''')


cursor.execute('''
INSERT INTO Courses (id, name, time_start, time_end)
VALUES (1, 'python', '21.07.21', '21.08.21');
''')
cursor.execute('''
INSERT INTO Courses (id, name, time_start, time_end)
VALUES (2, 'java', '13.07.21', '16.08.21');
''')


cursor.execute('''
INSERT INTO Students (id, name, surname, age, city)
VALUES (1, 'Max', 'Brooks', 24, 'Spb');
''')
cursor.execute('''
INSERT INTO Students (id, name, surname, age, city)
VALUES (2, 'John', 'Stones', 15, 'Spb');
''')
cursor.execute('''
INSERT INTO Students (id, name, surname, age, city)
VALUES (3, 'Andy', 'Wings', 45, 'Manhester');
''')
cursor.execute('''
INSERT INTO Students (id, name, surname, age, city)
VALUES (4, 'Kate', 'Brooks', 34, 'Spb');
''')


cursor.execute('''
INSERT INTO Student_courses (student_id, course_id)
VALUES (1, 1);
''')
cursor.execute('''
INSERT INTO Student_courses (student_id, course_id)
VALUES (2, 1);
''')
cursor.execute('''
INSERT INTO Student_courses (student_id, course_id)
VALUES (3, 1);
''')
cursor.execute('''
INSERT INTO Student_courses (student_id, course_id)
VALUES (4, 2);
''')


cursor.execute('''
SELECT * FROM Students
WHERE age > 30;
''')
students_over_30 = cursor.fetchall()
print("Студенты старше 30 лет:", students_over_30)


cursor.execute('''
SELECT Students.* FROM Students
JOIN Student_courses ON Students.id = Student_courses.student_id
JOIN Courses ON Courses.id = Student_courses.course_id
WHERE Courses.name = 'python';
''')
students_python = cursor.fetchall()
print("Студенты, проходящие курс по python:", students_python)


cursor.execute('''
SELECT Students.* FROM Students
JOIN Student_courses ON Students.id = Student_courses.student_id
JOIN Courses ON Courses.id = Student_courses.course_id
WHERE Courses.name = 'python' AND Students.city = 'Spb';
''')
students_python_spb = cursor.fetchall()
print("Студенты, проходящие курс по python и из Spb:", students_python_spb)


connection.commit()
connection.close()

print("База данных и таблицы успешно созданы и заполнены.")

